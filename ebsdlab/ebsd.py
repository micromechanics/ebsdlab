"""EBSD class: read, analyze and plot EBSD maps"""
import math
import time
from pathlib import Path
from typing import Any
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
from matplotlib import colormaps, colors
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from scipy.interpolate import griddata
from scipy.sparse import coo_array
from scipy.sparse.csgraph import connected_components
from scipy.spatial.transform import Rotation
from . import fileIO
from ._rotation import asBungeEulers, fromBungeEulers, multiply
from .symmetry import Symmetry, showOrSave


class EBSD:  # pylint: disable=too-many-public-methods
    """Class to allow for read EBSD data"""

    def __init__(self, fileName: str | Path, symmetry: str = '') -> None:
        """Initialize
        - Phases: phaseID 0 means not identified; phases are numbered from 1.
        - self.sym[k] is the symmetry of phase k, self.sym[0] is an empty Symmetry().

        Args:
           fileName: file name in the present directory
           symmetry: optional crystal symmetry name, e.g. "cubic".
               When supplied, it overrides the symmetries of all phases read from the file.
        """
        # initialize
        self.meta: dict[str, Any] = {}
        self.phi1: np.ndarray      = np.empty(0)
        self.phi: np.ndarray       = np.empty(0)
        self.phi2: np.ndarray      = np.empty(0)
        self.iq: np.ndarray        = np.empty(0)
        self.ci: np.ndarray        = np.empty(0)
        self.phaseID: np.ndarray   = np.empty(0)
        self.semSignal: np.ndarray = np.empty(0)
        self.fit: np.ndarray       = np.empty(0)
        # only some formats: Oxford band contrast/slope, bands, error, reliability index; OIM grain ID (or calcGrains)
        self.bc: np.ndarray        = np.empty(0)
        self.bs: np.ndarray        = np.empty(0)
        self.bands: np.ndarray     = np.empty(0)
        self.error: np.ndarray     = np.empty(0)
        self.ri: np.ndarray        = np.empty(0)
        self.grainID: np.ndarray   = np.empty(0)
        self.kam: np.ndarray       = np.empty(0)  # set by calcKAM
        self.cleaned: np.ndarray   = np.empty(0)  # set by grainDilation: points that were changed
        self.width     = 0.0
        self.height    = 0.0
        self.ratio     = 0.0
        self.stepSizeX = 0.0
        self.stepSizeY = 0.0  # row spacing
        # grid-coordinates
        self.grid        = 'SqrGrid'  # or 'HexGrid'
        self.nPoints     = 0
        self.nRows       = 0
        self.nColsOdd    = 0  # points in 1st, 3rd, ... row
        self.nColsEven   = 0  # points in 2nd, 4th, ... row
        self.x0, self.y0 = 0.0, 0.0
        self.xOffset     = 0.0  # x-shift of even rows: stepSizeX/2 for hexagonal grids
        startTime        = time.time()
        self.scanUnit    = 'um'
        self.sym: list[Symmetry] = []  # loaders add phase 1, 2, ...
        self._gridIndex: np.ndarray | None = None  # set by setGrid for unordered or partial maps

        # read input file header and parse it
        self.fileName = str(fileName)
        suffix        = Path(self.fileName).suffix.lower()
        if suffix in fileIO.LOADERS:
            fileIO.LOADERS[suffix](self)
        elif self.fileName.startswith('void'):
            print('Void mode', self.fileName[4:])
            fileIO.loadVoid(self, self.fileName[4:])
        else:
            raise ValueError('Unsupported file format. Supported formats: '+ ', '.join(set(fileIO.LOADERS)))

        # all files loaded
        if self._gridIndex is not None:  # place points onto the full grid; missing points: ci=-1, phaseID=0
            for name, values in list(vars(self).items()):
                if isinstance(values, np.ndarray) and values.shape == self._gridIndex.shape and name != '_gridIndex':
                    full = np.zeros(self.nPoints, dtype=values.dtype)
                    full[self._gridIndex] = values
                    setattr(self, name, full)
            missing = np.ones(self.nPoints, dtype=bool)
            missing[self._gridIndex] = False
            self.ci[missing] = -1
            print('   Missing points filled:', missing.sum())
            self._gridIndex = None
        if symmetry:
            self.sym = [Symmetry(symmetry)] * max(1, int(self.phaseID.max()))
        # phases without known symmetry are not identified
        self.sym = [Symmetry()] + self.sym + [Symmetry()] * (int(self.phaseID.max()) - len(self.sym))

        # output state
        print('   Read file with step size:', self.stepSizeX, self.stepSizeY, self.grid)
        print('   Optimal image pixel size:', int(self.width/self.stepSizeX))
        print('   Number of points:', self.nPoints)

        # convert into quaternions and only use that
        eulers = np.vstack((self.phi1, self.phi, self.phi2))
        self.quaternions = fromBungeEulers(eulers.T)
        fileIO.rotateToConventions(self, suffix)
        del self.phi1
        del self.phi
        del self.phi2

        # for plotting: determine image and imageSize once, use multiple times
        self.image: np.ndarray = np.empty(0)
        self.imageExtent       = (0.0, 0.0, 0.0, 0.0)  # of self.image: left, right, bottom, top in [um]
        self.mask              = self.ci > -1  # all are visible initially
        self.vMask             = np.ones(self.nPoints, dtype=bool)
        self.vMaskEvery        = 1  # preview: maps show only every k-th row and column
        print('   Duration init: ', int(np.round(time.time()-startTime)), 'sec')
        return


    ##
    # @name PLOT METHODS
    # @{
    def plot(self, vector: np.ndarray, widthPixel: int | None = None, vmax: float | None = None,
             vmin: float | None = None, interpolationType: str = 'nearest', cmap: Any = None,
             show: bool = True, cbar: bool = True) -> Any:
        """given a class-vector, plot the vector as an image

        Args:
           vector: vector to be plotted as a 2D image
           widthPixel: rescale to horizontal size of the image [default: optimal pixel width]
           vmax: rescale z-scale to maximal value
           vmin: rescale z-scale to minimal value
           interpolationType: interpolation type [default: "nearest" next-neighbor]
           cmap: colormap [default: Spectral with masked points in black]
           show: show the figure
           cbar: add a colorbar

        Returns:
           matplotlib figure
        """
        startTime = time.time()
        # create a special cmap palette with blacK as value for bad-numbers
        if cmap is None:
            cmap = colormaps['Spectral'].with_extremes(bad='k')
        z, self.imageExtent = self._image(vector, widthPixel, interpolationType)
        z = z.astype(float)
        mask, _ = self._image(~self.mask, widthPixel, interpolationType)
        # plot if/if-not the maximum and minimum are given
        fig, ax = plt.subplots()
        im = ax.imshow(np.ma.masked_where(mask, z), extent=self.imageExtent, cmap=cmap, vmax=vmax, vmin=vmin,
                       origin='upper')
        if cbar:
            fig.colorbar(im, ax=ax)
        print('   Plot with x and y axis in [um]')
        print('Duration plot: ', int(np.round(time.time()-startTime)), 'sec')
        if show:
            plt.show()
        z *= 255/np.max(z)
        self.image = z
        return fig


    def plotIPF(self, direction: str | int = 'ND', widthPixel: int | None = None, fileName: str = '',
                interpolationType: str = 'nearest', show: bool = True) -> Any:
        """
        plot Inverse Pole Figure (IPF)

        Args:
           direction: default.."ND", "RD", "TD"; a number is used as widthPixel with "ND"
           widthPixel: horizontal size of the image [default: optimal size based on data]
           fileName: save to file instead of showing
           interpolationType: interpolation type [default: "nearest"]
           show: show the figure if no fileName is given

        Returns:
           matplotlib figure
        """
        startTime = time.time()
        if direction == 'RD':
            axis = [1, 0, 0]
        elif direction == 'TD':
            axis = [0, 1, 0]
        elif direction == 'ND':
            axis = [0, 0, 1]
        else:  # if first argument specifies widthPixel
            widthPixel = int(direction)
            axis = [0, 0, 1]

        # colors only for the points shown in the image
        if interpolationType == 'nearest':
            shown = np.zeros(self.nPoints, dtype=bool)
            shown[self._image(np.arange(self.nPoints), widthPixel)[0]] = True
            shown = np.flatnonzero(shown)
        else:
            shown = np.flatnonzero(self.vMask)
        rgbs = np.zeros((3, self.nPoints), dtype=float)
        for phase, sym in enumerate(self.sym):
            points = shown[self.phaseID[shown] == phase]
            if not sym.lattice or not points.size:
                continue
            flags = np.zeros(len(points), dtype=bool)
            rgbsPhase = np.zeros((3, len(points)), dtype=float)
            # pole of the equivalent orientation q*s: (q*s)^-1 axis = s^-1 (q^-1 axis), so q^-1 axis is computed once
            poles = self.quaternions[points].inv().apply(axis).T
            for symmetry in sym.symmetryQuats():
                pole = symmetry.inv().as_matrix() @ poles[:, ~flags]
                remainingFlags, remainingRgbs = sym.inSST(pole, color=True, proper=False)
                if len(remainingRgbs.shape) == 2:
                    rgbsPhase[:, ~flags] = remainingRgbs
                    flags[~flags]        = remainingFlags
            rgbs[:, points] = rgbsPhase
        fig = self.plotRGB(rgbs, widthPixel, interpolationType)
        print('Duration plotIPF: ', int(np.round(time.time()-startTime)), 'sec')
        showOrSave(fileName, show)
        return fig


    def plotRGB(self, rgb: np.ndarray, widthPixel: int | None = None, interpolationType: str = 'nearest') -> Any:
        """given a RGB vector (same size as the other class vectors) plot the vector as an image

        Args:
           rgb: matrix [3, classVectorSize] to be plotted as a 2D image
           widthPixel: horizontal size of the image [default: optimal pixel width]
           interpolationType: interpolation type [default: "nearest"]

        Returns:
           matplotlib figure
        """
        # masked points are black
        rgb[:, ~self.mask] = 0
        image, self.imageExtent = self._image(rgb.T, widthPixel, interpolationType)
        self.image = (image*255).astype(np.uint8)
        fig, ax = plt.subplots()
        ax.imshow(self.image, extent=self.imageExtent, origin='upper')
        return fig


    def addSymbol(self, x: float, y: float, fileName: str = '', scale: float = 1., colorCube: str = 'black') -> None:
        """
        TODO: use version in ebsd_Orientation
        Add symbol of crystal orientation (symmetry and rotation) to IPF at given location

        Args:
           x: x-coordinate
           y: y-coordinate
           fileName: export to file
           scale: scale of symbol
           colorCube: color of symbol
        """
        # axes fill the figure exactly, so the rendered canvas is the image without margins
        fig = plt.figure(figsize=(6.4, 6.4/self.ratio), dpi=100)
        ax  = fig.add_axes((0, 0, 1, 1))
        ax.imshow(self.image, extent=self.imageExtent, origin='upper', aspect='auto')
        ax.axis('off')
        iClose = self.addUnitCellOverlay(ax, x, y, scale, colorCube)
        print('Euler angles at point:', np.round(asBungeEulers(self.quaternions[iClose], degrees=True), 1))
        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        self.image = np.asarray(canvas.buffer_rgba())[..., :3].copy()
        plt.close(fig)
        plt.imshow(self.image, extent=self.imageExtent, origin='upper')
        showOrSave(fileName)
        return


    def addUnitCellOverlay(self, ax: Any, x: float, y: float, scale: float = 1., colorCube: str = 'black') -> int:
        """
        Draw the nearest orientation's unit cell onto an existing IPF axis.

        This is the composable counterpart to :meth:`addSymbol`. It is useful
        for applications which manage the figure themselves, such as the GUI.

        Args:
           ax: matplotlib axis containing the IPF map
           x: x-coordinate in [um]
           y: y-coordinate in [um]
           scale: scale of unit cell
           colorCube: color of unit cell

        Returns:
           index of the point whose orientation is drawn
        """
        iClose      = int(self.nearestIndex(x, y))
        iQuaternion = self.quaternions[iClose]
        sym = self.sym[self.phaseID[iClose]]
        ax.autoscale(False)  # the overlay must not change the map limits
        if sym.lattice:
            # seen from above the sample (-Z): turned 180° about X, the edges toward the viewer have z > 0
            for start, end, lw in sym.unitCellSegments(Rotation.from_rotvec([np.pi, 0, 0])*iQuaternion, scale):
                ax.plot([x+start[0], x+end[0]], [y-start[1], y-end[1]], color=colorCube, lw=lw)
        return int(iClose)


    def addScaleBar(self, fileName: str = '', site: str = 'BL', barLength: float | None = None,
                    alpha: float = 0.5) -> Any:
        """
        Add scale-bar to image

        Args:
           fileName: if given, save to file
           site: where to put the scale bar: bottom-left "BL" (default), bottom-right "BR",
                 top-left "TL", top-right "TR"
           barLength: length of scale bar. It is calculated if not given
           alpha: transparency of scale bar background

        Returns:
           matplotlib figure
        """
        sites   = {'BL': 'lower left', 'BR': 'lower right', 'TL': 'upper left', 'TR': 'upper right'}
        fig, ax = plt.subplots()
        ax.imshow(self.image, extent=self.imageExtent, origin='upper')
        scaleBar = self.addScaleBarOverlay(ax, barLength, sites.get(site, 'lower left'))
        scaleBar.patch.set_alpha(alpha)
        ax.axis('off')
        showOrSave(fileName)
        return fig


    def addScaleBarOverlay(self, ax: Any, barLength: float | None = None, site: str = 'lower left') -> Any:
        """
        Add a scale bar to an existing map axis.

        Unlike :meth:`addScaleBar`, this preserves the supplied Matplotlib
        figure and is therefore suitable for interactive applications.

        Args:
           ax: Matplotlib axis containing an EBSD map in micrometres.
           barLength: scale-bar length in micrometres; calculated if omitted.
           site: Matplotlib anchored-artists location, e.g. ``"lower left"``.

        Returns:
           scale bar artist
        """
        if barLength is None:
            # visible area of the last image; the whole map before the first plot
            xLo, xHi, yHi, yLo = self.imageExtent if self.image.size else (0, self.width, self.height, 0)
            width, height = xHi-xLo, yHi-yLo
            digits    = int(math.log10(round(width/4.)))
            barLength = round(max(width, height) / 6., -digits)
        scaleBar = AnchoredSizeBar(ax.transData, barLength, str(barLength)+' '+'\u03BC'+'m', site, pad=0.5,
                                   color='black', frameon=True, size_vertical=barLength/15.,
                                  fontproperties=FontProperties(size=13.5))
        ax.add_artist(scaleBar)
        return scaleBar


    def plotPF(self, axis: Any = (1, 0, 0), points: bool = False,  # pylint: disable=too-many-locals
               fileName: str = '', color: str = '#1f77b4', alpha: float = 1.0, show: bool = True,
               density: int = 256, size: int = 2, vmin: float = 0.0, vmax: float | None = None,
               width: float = 5.0) -> Any:
        """plot pole figure, oriented as the map: X right, Y down, seen from above the sample (upper hemisphere -Z)

        The distribution is a pole density function in multiples of a random distribution (mrd): the poles are
        binned on an equal-area grid, smoothed with a von Mises-Fisher kernel on the sphere and then projected.

        Args:
          axis:    axis to plot: default: axis=1,0,0
          points:  plot individual points, or plot distribution [default]
          fileName: if given, save to file
          color:   plot color
          alpha:   alpha transparency
          show:    show figure [default], False for subsequent plotting
          density: distribution: size of the image in pixels
          size:    points: point size
          vmin:    distribution: lower values [mrd] are transparent
          vmax:    distribution: maximum of the color scale [mrd]; default: maximum of the data
          width:   distribution: half width of the kernel in degrees

        Returns:
          matplotlib figure
        """
        startTime = time.time()
        fig, ax = plt.subplots()
        maxColor = tuple(np.array(colors.hex2color(color))*0.5)
        axis = np.array(axis, dtype=float)/np.linalg.norm(axis)
        directions = []  # poles of all phases
        for phase, sym in enumerate(self.sym):
            if not sym.lattice:
                continue
            quaternions = self.quaternions[self.mask & self.vMask & (self.phaseID == phase)]
            for q in sym.equivalentQuaternions(Rotation.identity()):
                directions.append(quaternions.apply(q.apply(axis)))
        direction = np.concatenate(directions)
        if points:
            # upper hemisphere: toward the viewer above the sample; Y down in the plot
            direction = direction[direction[:, 2] < 0]
            x, y = direction[:, 0]/(1.-direction[:, 2]), -direction[:, 1]/(1.-direction[:, 2])
            ax.plot(x, y, '.', color=maxColor, markersize=size)
            ax.plot(np.cos(np.linspace(0, 2*np.pi, 100)),
                    np.sin(np.linspace(0, 2*np.pi, 100)), 'k-')
            ax.plot([-1, 1], [0, 0], 'k--')
            ax.plot([0, 0], [-1, 1], 'k--')
        else:
            # pole density function in multiples of a random distribution (mrd) on the upper hemisphere (-Z):
            # poles binned on an equal-area grid (Lambert projection, Y down), smoothed with the axial von Mises-Fisher
            # kernel exp(kappa (|cos| - 1)), where a direction and its opposite are the same pole, then projected
            kappa = np.log(2)/(1-np.cos(np.radians(width)))
            nGrid = int(np.clip(round(4*np.sqrt(2)/np.radians(width)), 48, 128))  # about 4 cells per half width
            direction = np.where(direction[:, 2:] > 0, -direction, direction)   # fold onto the upper hemisphere
            factor = np.sqrt(2/(1-direction[:, 2]))
            lambert = np.stack([-direction[:, 1]*factor, direction[:, 0]*factor])  # rows: plot y (Y down), columns: x
            cell = np.clip(((lambert + np.sqrt(2))/(2*np.sqrt(2))*nGrid).astype(int), 0, nGrid-1)
            counts = np.bincount(cell[0]*nGrid + cell[1], minlength=nGrid**2).astype(float)
            # cell centers on the sphere; centers outside the projection disk move onto its rim
            row, col = (np.indices((nGrid, nGrid)).reshape(2, -1) + 0.5)/nGrid*2*np.sqrt(2) - np.sqrt(2)
            rr = np.minimum(row**2 + col**2, 2.)
            centers = np.stack([col*np.sqrt(1-rr/4), -row*np.sqrt(1-rr/4), rr/2-1], axis=1)
            centers /= np.linalg.norm(centers, axis=1, keepdims=True)
            filled = counts > 0
            # ponytail: all pairs of cells, about 7 s for width 2°; a kd-tree with a cutoff angle if small widths matter
            pdf = np.zeros(nGrid**2)
            for start in range(0, nGrid**2, 1024):
                cosines = np.abs(centers[start:start+1024] @ centers[filled].T)
                pdf[start:start+1024] = np.exp(kappa*(cosines-1)) @ counts[filled]
            # uniform distribution: len(direction)/(2 pi) poles per steradian times the kernel integral (hemisphere)
            kernelIntegral = 2*np.pi*(1-np.exp(-kappa))/kappa
            pdf = (pdf * 2*np.pi/(len(direction)*kernelIntegral)).reshape(nGrid, nGrid)
            # stereographic image: pixel -> direction on the upper hemisphere -> Lambert grid of the density
            sx, sy = np.meshgrid(np.linspace(-1, 1, density), np.linspace(-1, 1, density))
            rr = sx**2 + sy**2
            scale = np.sqrt(2/(1+(1-rr)/(1+rr)))*2/(1+rr)  # Lambert coordinates are the stereographic ones times this
            rows, cols = ((c*scale + np.sqrt(2))/(2*np.sqrt(2))*nGrid - 0.5 for c in (sy, sx))
            img = ndi.map_coordinates(pdf, [rows, cols], order=1, mode='nearest')
            img[(rr > 1) | (img < vmin)] = np.nan
            cmap = colors.LinearSegmentedColormap.from_list('my', [(1, 1, 1), maxColor])
            image = ax.imshow(img, cmap=cmap, alpha=alpha, vmin=0.0, vmax=vmax, origin='lower', extent=(-1, 1, -1, 1))
            fig.colorbar(image, ax=ax, label='mrd', shrink=0.8)
            ax.plot(np.cos(np.linspace(0, 2*np.pi, 100)),
                    np.sin(np.linspace(0, 2*np.pi, 100)), 'k-', lw=2)
            ax.plot([0, 0], [-1, 1], 'k--', lw=1)
            ax.plot([-1, 1], [0, 0], 'k--', lw=1)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim((-1, 1))
        ax.set_ylim((-1, 1))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        print('Duration plotPF: ', int(np.round(time.time()-startTime)), 'sec')
        showOrSave(fileName, show)
        return fig

    # @}


    ##
    # @name Mask and Path routines
    # Masked areas are plotted in black. Initially no point is part of the mask;
    # all points are false.
    # @{
    def maskCI(self, ci: float) -> None:
        """masked all points off, which have a CI less than: good points=False, bad points=True

        Args:
           ci: critical CI
        """
        self.mask = self.ci > ci
        return


    def maskReset(self) -> None:
        """reset mask"""
        self.mask = self.ci > -1
        return


    def removePointsOutsideMask(self) -> None:
        """Set points outside the mask invalid so they are read as nonexistent by OIM after export"""
        self.ci[~self.mask] = -1.0
        self.fit[~self.mask] = 180.0
        return


    def setVMask(self, every: int = 1) -> None:
        """fast preview, does not influence results in any way: maps show only every k-th row and column
        (at any image size), pole figures use every k-th point.

        Args:
           every: use only every k-th point. Improves plotting speed. every=1 resets.
        """
        self.vMask[:] = False
        self.vMask[::every] = True
        self.vMaskEvery = every
        return


    def cropVMask(self, xmin: float | None = None, ymin: float | None = None,
                  xmax: float | None = None, ymax: float | None = None) -> None:
        """crop visible area

        Args:
           xmin: minimum x-coordinate
           ymin: minimum y-coordinate
           xmax: maximum x-coordinate
           ymax: maximum y-coordinate
        """
        xmin, ymin = xmin or 0, ymin or 0
        x, y = self.xy()
        xmax, ymax = xmax or np.max(x), ymax or np.max(y)
        self.vMask = self.vMask & (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
        return


    def xy(self, idx: Any = None) -> tuple[np.ndarray, np.ndarray]:
        """x- and y-coordinates of points, computed from the grid parameters

        Args:
           idx: point indices [if None: all points]

        Returns:
           x, y
        """
        i = np.arange(self.nPoints) if idx is None else np.asarray(idx)
        pair, j = np.divmod(i, self.nColsOdd + self.nColsEven)
        even = j >= self.nColsOdd  # 2nd row of the pair of rows
        col = j - even*self.nColsOdd
        return self.x0 + col*self.stepSizeX + even*self.xOffset, self.y0 + (2*pair+even)*self.stepSizeY


    def calcKAM(self) -> None:
        """calculate Kernel Average Misorientation in DEGREES (because user focused) from the nearest neighbors"""
        startTime = time.time()
        _, angles = self._neighborMisorientations()
        self.kam = np.degrees(np.nanmean(angles, axis=1))
        self.kam[self.ci == -1.0] = np.nan
        print('Duration KAM evaluation: ', int(np.round(time.time()-startTime)), 'sec')
        return


    def calcGrains(self, tolerance: float = 5., minSize: int = 6, minNRows: int = 2) -> None:
        """identify grains: neighbors of one phase that are indexed and misoriented less than tolerance are linked;
        grains are the connected regions. grainID is 1, 2, ...; 0 = unindexed or grain too small.

        Args:
           tolerance: maximum misorientation of neighbors in a grain in degrees
           minSize: minimum number of points of a grain
           minNRows: minimum number of rows and of columns a grain spans
        """
        startTime = time.time()
        neighbors, angles = self._neighborMisorientations()
        linked = angles < np.radians(tolerance)  # nan: not linked
        points = np.broadcast_to(np.arange(self.nPoints)[:, None], neighbors.shape)
        graph  = coo_array((np.ones(linked.sum()), (points[linked], neighbors[linked])),
                           shape=(self.nPoints, self.nPoints))
        _, labels = connected_components(graph, directed=False)
        x, y = self.xy()
        rows = np.rint((y-self.y0)/self.stepSizeY).astype(int)
        cols = np.floor((x-self.x0)/self.stepSizeX + 0.01).astype(int)  # hex: the shifted rows share a column
        keep = np.bincount(labels) >= minSize
        for coordinate in (rows, cols):  # number of distinct rows/columns of each region
            period = coordinate.max()+1
            keep &= np.bincount(np.unique(labels*period + coordinate)//period, minlength=len(keep)) >= minNRows
        indexed = np.array([bool(sym.lattice) for sym in self.sym])[self.phaseID] & (self.ci != -1.0)
        keep[labels[~indexed]] = False  # unindexed points are single regions
        self.grainID = (np.cumsum(keep)*keep).astype(np.uint32)[labels]
        print('   Number of grains:', int(keep.sum()))
        print('Duration grain identification: ', int(np.round(time.time()-startTime)), 'sec')
        return


    def grainDilation(self) -> None:
        """fill points without grain (grainID 0) from the neighbors, in place (reload to undo); runs calcGrains with
        default arguments if it has not run. A point joins the grain of the majority of its neighbors that have a grain
        and takes phase, orientation, CI and mask of the neighbor with the highest CI among them; this also decides
        ties.
        Repeat until nothing changes. Changed points are True in self.cleaned.
        """
        startTime = time.time()
        if not self.grainID.size:
            self.calcGrains()
        neighbors = self.neighbors()
        # missing neighbors (-10): the point itself, which has no grain and therefore no vote
        neighbors = np.where(neighbors >= 0, neighbors, np.arange(self.nPoints)[:, None])
        self.cleaned = self.grainID == 0
        while True:
            grains  = self.grainID[neighbors]
            points  = np.flatnonzero((self.grainID == 0) & (grains > 0).any(axis=1))
            if not points.size:
                break
            grains  = grains[points]
            votes   = ((grains[:, :, None] == grains[:, None, :]) & (grains[:, None, :] > 0)).sum(axis=2)
            ciVotes = np.where((grains > 0) & (votes == votes.max(axis=1, keepdims=True)),
                               self.ci[neighbors[points]], -np.inf)
            source  = neighbors[points, np.argmax(ciVotes, axis=1)]
            self.grainID[points]     = self.grainID[source]
            self.phaseID[points]     = self.phaseID[source]
            self.quaternions[points] = self.quaternions[source]
            self.ci[points]          = self.ci[source]
            self.mask[points]        = self.mask[source]
        self.cleaned &= self.grainID > 0
        print('   Number of points changed:', int(self.cleaned.sum()))
        print('Duration grain dilation: ', int(np.round(time.time()-startTime)), 'sec')
        return


    def _neighborMisorientations(self) -> tuple[np.ndarray, np.ndarray]:
        """misorientation of each point to its nearest neighbors

        Returns:
           neighbors (see neighbors()), angles in radians; nan for invalid, unindexed or other-phase neighbors
        """
        neighbors = self.neighbors()
        # for indexing, missing neighbors (-10) are the point itself; their angles are nan below
        safe   = np.where(neighbors >= 0, neighbors, np.arange(len(neighbors))[:, None])
        angles = np.full(neighbors.shape, np.nan)
        for phase, sym in enumerate(self.sym):
            points = np.flatnonzero(self.phaseID == phase)
            if not sym.lattice or not points.size:
                continue
            # scalar part of misQ*s for all symmetries s: dot product with (-s_x, -s_y, -s_z, s_w)
            symQ = sym.symmetryQuats().as_quat() * [-1, -1, -1, 1]
            for iNeighbor in range(neighbors.shape[1]):
                misQ = multiply(self.quaternions[points].inv(), self.quaternions[safe[points, iNeighbor]])
                scalar = np.abs(misQ.as_quat() @ symQ.T).max(axis=1)
                angles[points, iNeighbor] = 2*np.arccos(np.clip(scalar, 0., 1.))
        angles[(neighbors < 0) | (self.phaseID[safe] != self.phaseID[:, None])] = np.nan
        angles[self.ci[safe] == -1.0] = np.nan
        angles[self.ci == -1.0] = np.nan
        return neighbors, angles


    def neighbors(self, idx: int | None = None) -> np.ndarray:
        """identify the nearest neighboring indexes

        Args:
           idx: index to find [if None: calculate all]

        Returns:
         array of neighbors; invalid points have a value=-10
        """
        i        = np.arange(self.nPoints) if idx is None else np.atleast_1d(idx)
        period   = self.nColsOdd + self.nColsEven
        pair, j  = np.divmod(i, period)
        even     = (j >= self.nColsOdd).astype(int)  # 2nd row of the pair of rows
        row, col = 2*pair + even, j - even*self.nColsOdd
        if self.grid == 'HexGrid':
            # first of the two touching columns in the rows above and below; depends on shift direction
            shift = int(self.xOffset > 0)
            lower = col + np.where(even, shift-1, -shift)
            nRow  = row[:, None] + np.array([-1, -1, 0, 0, 1, 1])
            nCol  = np.stack([lower, lower+1, col-1, col+1, lower, lower+1], axis=1)
        else:
            nRow  = row[:, None] + np.array([-1, 0, 0, 1])
            nCol  = np.stack([col, col-1, col+1, col], axis=1)
        rowLength = np.where(nRow % 2, self.nColsEven, self.nColsOdd)
        neighbors = (nRow//2)*period + (nRow % 2)*self.nColsOdd + nCol
        valid     = (nRow >= 0) & (nRow < self.nRows) & (nCol >= 0) & (nCol < rowLength) & (neighbors < self.nPoints)
        neighbors[~valid] = -10
        if idx is not None:
            return neighbors[0]
        return neighbors
    # @}


    def writeANG(self, fileName: str) -> None:
        """WRAPPER: write body of ang file

        Args:
           fileName: file name
        """
        fileIO.writeANG(self, fileName)


    @property
    def x(self) -> np.ndarray:
        """x-coordinates of all points

        Returns:
           x-coordinates in [um]
        """
        return self.xy()[0]


    @property
    def y(self) -> np.ndarray:
        """y-coordinates of all points

        Returns:
           y-coordinates in [um]
        """
        return self.xy()[1]


    def setGrid(self, x: np.ndarray, y: np.ndarray) -> None:
        """Derive grid parameters from the coordinates; the coordinates themselves are not stored.
        If the points are unordered, partial (e.g. partition exports) or a single row, self._gridIndex is the grid
        position of each point; __init__ then places all per-point data onto the full grid.

        Args:
           x: x-coordinates of all points
           y: y-coordinates of all points
        """
        self._gridIndex = None
        rowStarts = np.flatnonzero(np.diff(x) < 0) + 1
        if len(rowStarts):  # complete row-major grid
            self.nPoints     = len(x)
            self.nRows       = len(rowStarts) + 1
            self.nColsOdd    = int(rowStarts[0])
            self.nColsEven   = int(rowStarts[1] if len(rowStarts) > 1 else len(x)) - self.nColsOdd
            self.x0, self.y0 = float(x[0]), float(y[0])
            self.stepSizeX   = float(x[1] - x[0])
            self.stepSizeY   = float(y[self.nColsOdd] - y[0])
            self.xOffset     = float(x[self.nColsOdd] - x[0])
            self.grid        = 'HexGrid' if abs(self.xOffset) > self.stepSizeX/10 else 'SqrGrid'
            xGrid, yGrid     = self.xy()
            if np.allclose(xGrid, x, atol=self.stepSizeX/10) and np.allclose(yGrid, y, atol=self.stepSizeY/10):
                return

        def smallestStep(gaps: np.ndarray) -> float:
            gaps = gaps[gaps > 1e-3*gaps.max()] if len(gaps) else gaps  # ignore rounding noise
            return float(gaps.min()) if len(gaps) else 0.0
        self.y0        = float(y.min())
        self.stepSizeY = smallestStep(np.diff(np.unique(y)))
        row            = np.rint((y-self.y0)/self.stepSizeY).astype(int) if self.stepSizeY else np.zeros(len(y), int)
        order          = np.lexsort((x, row))
        self.stepSizeX = smallestStep(np.diff(x[order])[np.diff(row[order]) == 0])
        self.stepSizeX = self.stepSizeX or self.stepSizeY  # single column
        self.stepSizeY = self.stepSizeY or self.stepSizeX  # single row
        if not self.stepSizeX:
            raise ValueError('EBSD data must contain more than one point.')
        even         = row % 2 == 1  # 2nd row of the pair of rows
        self.xOffset = 0.0
        if even.any() and abs((x[even].min()-x[~even].min())/self.stepSizeX % 1 - 0.5) < 0.1:
            self.xOffset = self.stepSizeX/2
        self.grid      = 'HexGrid' if self.xOffset else 'SqrGrid'
        self.x0        = float(min(x[~even].min(), x[even].min()-self.xOffset if even.any() else np.inf))
        col            = np.rint((x - self.x0 - even*self.xOffset)/self.stepSizeX).astype(int)
        self.nRows     = int(row.max()) + 1
        if self.grid == 'HexGrid':
            self.nColsOdd, self.nColsEven = int(col[~even].max()) + 1, int(col[even].max()) + 1
        else:
            self.nColsOdd  = int(col.max()) + 1
            self.nColsEven = self.nColsOdd if self.nRows > 1 else 0
        period          = self.nColsOdd + self.nColsEven
        self.nPoints    = (self.nRows//2)*period + (self.nRows % 2)*self.nColsOdd
        self._gridIndex = (row//2)*period + even*self.nColsOdd + col
        return


    def _image(self, values: np.ndarray, widthPixel: int | None = None,
               interpolationType: str = 'nearest') -> tuple[np.ndarray, tuple[float, float, float, float]]:
        """Sample per-point values onto the pixels of an image of the visible area

        Args:
           values: one value (or one row of values) per point
           widthPixel: horizontal size of the image [default: one pixel per step; half a step for hexagonal grids]
           interpolationType: "nearest" uses the grid directly; others interpolate with scipy's griddata

        Returns:
           image of shape (height, width, ...), extent for imshow in [um]
        """
        xVisible, yVisible = self.xy(np.flatnonzero(self.vMask))
        xLo, xHi           = xVisible.min()-self.stepSizeX/2, xVisible.max()+self.stepSizeX/2
        yLo, yHi           = yVisible.min()-self.stepSizeY/2, yVisible.max()+self.stepSizeY/2
        self.ratio         = (xHi-xLo) / (yHi-yLo)
        if widthPixel is None:
            pitch          = self.stepSizeX/2 if self.grid == 'HexGrid' else self.stepSizeX
            widthPixel     = int(round((xHi-xLo)/pitch))
            heightPixel    = int(round((yHi-yLo)/self.stepSizeY))
        else:
            heightPixel    = int(round(widthPixel/self.ratio))
        xAxis              = xLo + (np.arange(widthPixel)+0.5)*(xHi-xLo)/widthPixel
        yAxis              = yLo + (np.arange(heightPixel)+0.5)*(yHi-yLo)/heightPixel
        extent             = (xLo, xHi, yHi, yLo)
        if interpolationType != 'nearest':
            x, y = np.meshgrid(xAxis, yAxis)
            return griddata(np.vstack((xVisible, yVisible)).T, values[self.vMask], (x, y), interpolationType), extent
        if self.vMaskEvery > 1:  # preview: pixels sample only every k-th row and column
            xBlock, yBlock = self.stepSizeX*self.vMaskEvery, self.stepSizeY*self.vMaskEvery
            xAxis          = xLo + (np.floor((xAxis-xLo)/xBlock)+0.5)*xBlock
            yAxis          = yLo + (np.floor((yAxis-yLo)/yBlock)+0.5)*yBlock
        return values[self.nearestIndex(xAxis[None, :], yAxis[:, None])], extent


    def nearestIndex(self, x: Any, y: Any) -> np.ndarray:
        """Index of the grid point closest to x, y: row from y, then column within that row from x

        Args:
           x: x-coordinate(s) in [um]; broadcast with y
           y: y-coordinate(s) in [um]

        Returns:
           point indices, of the broadcast shape of x and y
        """
        row  = np.clip(np.rint((np.asarray(y)-self.y0)/self.stepSizeY).astype(int), 0, self.nRows-1)
        even = row % 2
        col  = np.rint((np.asarray(x) - self.x0 - even*self.xOffset)/self.stepSizeX).astype(int)
        col  = np.clip(col, 0, np.where(even, self.nColsEven, self.nColsOdd)-1)
        idx  = (row//2)*(self.nColsOdd+self.nColsEven) + even*self.nColsOdd + col
        return np.minimum(idx, self.nPoints-1)
