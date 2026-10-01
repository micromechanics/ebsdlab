##
# @file
# @brief Class to allow for read EBSD data
#
import math
import os
import time
from pathlib import Path
from typing import Any
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
from matplotlib import colors
from matplotlib.backends.backend_agg import FigureCanvasAgg
from scipy.interpolate import griddata
from scipy.spatial.transform import Rotation
from ._rotation import asBungeEulers
from .orientation import Orientation
from .symmetry import Symmetry

SUPPORTED_SUFFIXES = {'.ang', '.osc', '.txt', '.crc'}
TSL_SYMMETRIES = {43: 'cubic', 62: 'hexagonal', 42: 'tetragonal', 22: 'orthorhombic', 32: 'trigonal',
                  2: 'monoclinic', 1: 'triclinic', 'm-3m': 'cubic'}
OXFORD_LAUE_GROUPS = {11: 'cubic', 9: 'hexagonal', 5: 'tetragonal', 3: 'orthorhombic', 7: 'trigonal',
                      2: 'monoclinic', 1: 'triclinic'}


class EBSD:
    """Class to allow for read EBSD data.

    Uses quaternions and symmetries but not orientations
    """

    def __init__(self, fileName: str | Path, symmetry: str = '') -> None:
        """
        read input file <br>
        initialize things<br>
        .ang or .osc file format
        - Header: ASCII information starting with #

        Args:
           fileName: file name in the present directory
           symmetry: optional crystal symmetry name, e.g. "cubic".
               When supplied, it overrides the symmetries of all phases read from the file.

        Phases: phaseID 0 means not identified; phases are numbered from 1. self.sym[k] is the
        symmetry of phase k, self.sym[0] is an empty Symmetry().
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
        self.width     = 0.0
        self.height    = 0.0
        self.ratio     = 0.0
        self.stepSizeX = 0.0
        self.stepSizeY = 0.0  # row spacing
        # grid: coordinates are computed from these instead of being stored, see xy()
        self.grid = 'SqrGrid'  # or 'HexGrid'
        self.nPoints = 0
        self.nRows = 0
        self.nColsOdd = 0   # points in 1st, 3rd, ... row
        self.nColsEven = 0  # points in 2nd, 4th, ... row
        self.x0, self.y0 = 0.0, 0.0
        self.xOffset = 0.0  # x-shift of even rows: stepSizeX/2 for hexagonal grids
        startTime = time.time()
        self.scanUnit = 'um'
        self.sym: list[Symmetry] = []  # loaders add phase 1, 2, ...

        # read input file header and parse it
        self.fileName = str(fileName)
        suffix = Path(self.fileName).suffix.lower()
        if suffix == '.ang':
            self.loadANG()
        elif suffix == '.osc':
            self.loadOSC()
        elif suffix == '.txt':
            self.loadTXT()
        elif suffix == '.crc':
            self.loadCRC()
        elif self.fileName.startswith('void'):
            print('Void mode', self.fileName[4:])
            self.loadVoid(self.fileName[4:])
        else:
            raise ValueError(
                'Unsupported EBSD file format. Supported formats are '
                + ', '.join(sorted(SUPPORTED_SUFFIXES)) + '.')

        if symmetry:
            self.sym = [Symmetry(symmetry)] * max(1, int(self.phaseID.max()))
        # phases without known symmetry are not identified
        self.sym = [Symmetry()] + self.sym + [Symmetry()] * (int(self.phaseID.max()) - len(self.sym))

        print('   Read file with step size:', self.stepSizeX, self.stepSizeY, self.grid)
        print('   Optimal image pixel size:', int(self.width/self.stepSizeX))
        print('   Number of points:', self.nPoints)

        # convert into quaternions and only use that
        eulers = np.vstack((self.phi1, self.phi, self.phi2))
        self.quaternions = Rotation.from_euler('ZXZ', eulers.T)
        del self.phi1
        del self.phi
        del self.phi2

        # for plotting: determine image and imageSize once, use multiple times
        self.image: np.ndarray = np.empty(0)
        self.imageExtent = (0.0, 0.0, 0.0, 0.0)  # of self.image: left, right, bottom, top in [um]
        self.mask = self.ci > -1  # all are visible initially
        self.vMask = np.ones(self.nPoints, dtype=bool)
        self.vMaskEvery = 1  # preview: maps show only every k-th row and column
        print('   Duration init: ', int(np.round(time.time()-startTime)), 'sec')
        return


    ##
    # @name PLOT METHODS
    # @{
    def plot(self, vector: np.ndarray, widthPixel: int | None = None, vmax: float | None = None,
             vmin: float | None = None, interpolationType: str = 'nearest', cmap: Any = None,
             show: bool = True, cbar: bool = True) -> Any:
        """
        given a class-vector, plot the vector as an image<br>
        the x and y are given by the grid

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
            cmap = cm.Spectral
            cmap.set_bad('k', 1.0)
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
            shown = np.unique(self._image(np.arange(self.nPoints), widthPixel)[0])
        else:
            shown = np.flatnonzero(self.vMask)
        rgbs = np.zeros((3, self.nPoints), dtype=float)
        for phase, sym in enumerate(self.sym):
            points = shown[self.phaseID[shown] == phase]
            if not sym.lattice or not len(points):
                continue
            flags = np.zeros(len(points), dtype=bool)
            rgbsPhase = np.zeros((3, len(points)), dtype=float)
            equivQuaternions = sym.equivalentQuaternions(self.quaternions[points])
            for equivQuaternion in equivQuaternions:
                pole = equivQuaternion.inv().apply(axis)
                remainingFlags, remainingRgbs = sym.inSST(
                    pole[~flags].T, color=True, proper=False)
                if len(remainingRgbs.shape) == 2:
                    rgbsPhase[:, ~flags] = remainingRgbs
                    flags[~flags] = remainingFlags
            rgbs[:, points] = rgbsPhase
        fig = self.plotRGB(rgbs, widthPixel, interpolationType)
        print('Duration plotIPF: ', int(np.round(time.time()-startTime)), 'sec')
        if not fileName and show:
            plt.show()
        elif fileName:
            plt.savefig(fileName, dpi=150, bbox_inches='tight')
            plt.close()
        return fig


    def plotRGB(self, rgb: np.ndarray, widthPixel: int | None = None, interpolationType: str = 'nearest') -> Any:
        """
        given a RGB vector (same size as the other class vectors)
        plot the vector as an image<br>
        the x and y are given by the grid
        USED INTERNALLY

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
        ax = fig.add_axes((0, 0, 1, 1))
        ax.imshow(self.image, extent=self.imageExtent, origin='upper', aspect='auto')
        ax.axis('off')
        iClose = self.addUnitCellOverlay(ax, x, y, scale, colorCube)
        print('Euler angles at point:',
              np.round(asBungeEulers(self.quaternions[iClose], degrees=True), 1))
        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        self.image = np.asarray(canvas.buffer_rgba())[..., :3].copy()
        plt.close(fig)
        plt.imshow(self.image, extent=self.imageExtent, origin='upper')
        if not fileName:
            plt.show()
        else:
            plt.savefig(fileName, dpi=150, bbox_inches='tight')
            plt.close()
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
        iClose = int(self._nearestIndex(x, y))
        iQuaternion = self.quaternions[iClose]
        sym = self.sym[self.phaseID[iClose]]
        if sym.lattice:
            for line in sym.unitCell():
                start = iQuaternion.apply(np.array(line[:3], dtype=float)*scale)
                end = iQuaternion.apply(np.array(line[3:], dtype=float)*scale)
                # OIM coordinate system and ``imshow(origin='upper')``.
                start = np.array([-start[1], -start[0], start[2]])
                end = np.array([-end[1], -end[0], end[2]])
                if start[2] < 0 and end[2] < 0:
                    segments = [(start, end, 0.2)]
                elif start[2] > 0 and end[2] > 0:
                    segments = [(start, end, 2)]
                else:
                    delta = end-start
                    mid = start+(-start[2]/delta[2])*delta
                    if start[2] > 0:
                        segments = [(start, mid, 2), (mid, end, 0.2)]
                    else:
                        segments = [(start, mid, 0.2), (mid, end, 2)]
                for first, last, lw in segments:
                    ax.plot([first[0]+x, last[0]+x], [first[1]+y, last[1]+y], color=colorCube, lw=lw)
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
        sites = {'BL': 'lower left', 'BR': 'lower right', 'TL': 'upper left', 'TR': 'upper right'}
        fig, ax = plt.subplots()
        ax.imshow(self.image, extent=self.imageExtent, origin='upper')
        scaleBar = self.addScaleBarOverlay(ax, barLength, sites.get(site, 'lower left'))
        scaleBar.patch.set_alpha(alpha)
        ax.axis('off')
        if not fileName:
            plt.show()
        else:
            plt.savefig(fileName, dpi=150, bbox_inches='tight')
            plt.close()
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
        from matplotlib.font_manager import FontProperties
        from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

        if barLength is None:
            # visible area of the last image; the whole map before the first plot
            xLo, xHi, yHi, yLo = self.imageExtent if self.image.size else (0, self.width, self.height, 0)
            width, height = xHi-xLo, yHi-yLo
            digits = int(math.log10(round(width/4.)))
            barLength = round(max(width, height) / 6., -digits)
        scaleBar = AnchoredSizeBar(ax.transData, barLength,
                                  str(barLength)+' '+'\u03BC'+'m', site,
                                  pad=0.5, color='black', frameon=True,
                                  size_vertical=barLength/15.,
                                  fontproperties=FontProperties(size=13.5))
        ax.add_artist(scaleBar)
        return scaleBar


    def plotPF(self, axis: Any = (1, 0, 0), points: bool = False, fileName: str = '',
               color: str = '#1f77b4', alpha: float = 1.0, show: bool = True, density: int = 256, size: int = 2,
               proj2D: str = 'up-left', vmin: float = 0.0, vmax: float = 1.0) -> Any:
        """
        plot pole figure

        Projection onto 2D: cooradinate systems are given as xDirection-yDirection (z follows)
        - down-right: [default in text books, mTex] RD = x = down; TD = y = right; ND = z = outOfPlane
        - up-left: [default in OIM and here] RD = x = up; TD = y = left; ND = z = outOfPlane

        Args:
          axis:    axis to plot: default: axis=1,0,0
          points:  plot individual points [default], or plot distribution
          fileName: if given, save to file
          color:   plot color
          alpha:   alpha transparency
          show:    show figure [default], False for subsequent plotting
          density: how many points to plot on the distribution
          size:    points: point size; distribution: amount of smoothing: higher more smoothing
          proj2D:  orientation of 2D projection: [down-right, up-left, None]
          vmin:    minimum value plotted, used as cut-off for transparency
          vmax:    max. used in color coding, allows to focus on minor texture

        Returns:
          matplotlib figure; None for an unknown proj2D
        """
        startTime = time.time()
        fig, ax = plt.subplots()
        maxColor = tuple(np.array(colors.hex2color(color))*0.5)
        for phase, sym in enumerate(self.sym):
            if not sym.lattice:
                continue
            oHelp = Orientation(eulers=np.array([0., 0., 0.]), symmetry=sym.lattice)
            axis = np.array(axis, dtype=float)
            axis /= np.linalg.norm(axis)
            mask = self.mask & self.vMask & (self.phaseID == phase)
            xs, ys = [], []
            for q in oHelp.symmetry.equivalentQuaternions(oHelp.quaternion):
                conjAxis = q.apply(axis)
                direction = self.quaternions.apply(conjAxis)
                direction = direction[mask]  # filter mask
                # filter upward dome
                direction = direction[direction[:, 2] > 0]
                direction[:, 0] /= direction[:, 2]+1.
                direction[:, 1] /= direction[:, 2]+1.
                xs.append(direction[:, 0])
                ys.append(direction[:, 1])
        x, y = np.concatenate(xs), np.concatenate(ys)
        if points:
            if proj2D == 'down-right':
                ax.plot(-x, y, '.', color=maxColor,
                        markersize=size)  # markersize=0.05
            elif proj2D == 'up-left':
                ax.plot(-y, x, '.', color=maxColor,
                        markersize=size)  # markersize=0.05
            else:
                return
            ax.plot(np.cos(np.linspace(0, 2*np.pi, 100)),
                    np.sin(np.linspace(0, 2*np.pi, 100)), 'k-')
            ax.plot([-1, 1], [0, 0], 'k--')
            ax.plot([0, 0], [-1, 1], 'k--')
        else:
            cmap = colors.LinearSegmentedColormap.from_list(
                'my', [(1, 1, 1), maxColor])
            center = (density - 1)/2
            imgDim = density+2*size
            img = np.zeros((imgDim, imgDim))
            x, y = np.nan_to_num(x), np.nan_to_num(y)
            if proj2D == 'down-right':
                zippedList = list(zip(-x, y))
            elif proj2D == 'up-left':
                zippedList = list(zip(-y, x))
            else:
                return
            for xCoordinate, yCoordinate in zippedList:
                ix = int((xCoordinate - -1.) * center) + size
                iy = int((yCoordinate - -1.) * center) + size
                if 0 <= ix < imgDim and 0 <= iy < imgDim:
                    img[iy][ix] += 1
            img = ndi.gaussian_filter(
                img, (size, size))  # gaussian convolution
            img /= np.max(img)                               # normalize
            # filter out low values to make transparent
            img[img < vmin] = np.nan
            ax.imshow(img, cmap=cmap, alpha=alpha,
                      vmin=0.0, vmax=vmax, origin='lower',
                      extent=(-1, 1, -1, 1))
            ax.plot(np.cos(np.linspace(0, 2*np.pi, 100)),
                    np.sin(np.linspace(0, 2*np.pi, 100)), 'k-', lw=2)
            ax.plot([0, 0], [-1, 1], 'k--', lw=1)
            ax.plot([-1, 1], [0, 0], 'k--', lw=1)
            # plt.colorbar()
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim((-1, 1))
        ax.set_ylim((-1, 1))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        print('Duration plotPF: ', int(np.round(time.time()-startTime)), 'sec')
        if not fileName and show:
            plt.show()
        elif fileName:
            plt.savefig(fileName, dpi=150, bbox_inches='tight')
            plt.clf()
            plt.cla()
        return fig

    # @}


    ##
    # @name INPUT METHODS
    # @{
    def loadANG(self, fileName: str = '') -> None:
        """
        Load .ang file: filename saved in self. No need to use it

        Args:
           fileName: file to read [default: self.fileName]
        """
        if fileName:
            self.fileName = fileName
        print('Load .ang file: ', self.fileName)
        keys = ['MaterialName', 'LatticeConstants',
                'WorkingDistance', 'SEMVoltage', 'GRID:', 'Symmetry']
        fileHandle = open(self.fileName)
        keyValues: list[Any] = [''] * len(keys)  # actual values
        symmetries = []
        for line in fileHandle:
            if line[0:10] == '# OPERATOR':
                break
            for key in keys:
                searchTerm = '# '+key
                if searchTerm == line[0:len(searchTerm)]:
                    index = keys.index(key)
                    value: Any = line.rstrip().split()[2:]
                    if len(value) == 1:
                        value = value[0]
                        try:
                            value = float(value)
                        except ValueError:
                            pass
                    keyValues[index] = value
                    if key == 'Symmetry':
                        symmetries.append(value)
                    break
        self.meta = dict(list(zip(keys, keyValues)))
        self.sym = [Symmetry(TSL_SYMMETRIES.get(i, '')) for i in symmetries]
        # read data: print "Reading file, this can take a bit..."
        data = np.loadtxt(fileHandle)
        self.phi1 = data[:, 0].astype(float)
        self.phi = data[:, 1].astype(float)
        self.phi2 = data[:, 2].astype(float)
        self.iq = data[:, 5].astype(float)
        self.ci = data[:, 6].astype(float)
        self.phaseID = data[:, 7].astype(np.uint8)
        self.phaseID += not self.phaseID.any()  # single-phase EDAX files use 0
        self.semSignal = data[:, 8].astype(np.uint8)
        self.fit = data[:, 9].astype(float)
        self.width  = max(data[:, 3])
        self.height = max(data[:, 4])
        self.ratio = self.width/self.height
        self._setGrid(data[:, 3], data[:, 4])
        fileHandle.close()
        del data
        return


    def loadTXT(self, fileName: str = '', update: bool = False) -> None:
        """
        read txt file and possibly update data. Warning, this resets the mask to the one of the file

        Update makes more sense if you have original data and update it with some partial information,
        since the partial information is incomplete (stepSize, width and height are in many cases wrong).
        Update will keep the old data and overwrite the new. THIS IS SLOWER THAN CREATING NEW

        Args:
          fileName: fileName to load (partition data from OIM)
          update: update data or read new (read-new: default)
        """
        print('TODO: Symmetry has to be read and used')
        startTime = time.time()
        print('Load .txt file:', fileName)
        if not fileName:
            fileName = self.fileName
        fileHandle = open(fileName)
        foundKeys: dict[str, int] = {}
        for line in fileHandle:
            if line[0] != '#':
                break
            parts = line.split()
            if len(parts) < 2:
                continue
            if parts[1] == 'Header:':
                print('   Header: ', parts[2])
                continue
            if 'Column' in parts:
                foundKeys[parts[3]] = int(parts[2].split(':')[0].split('-')[0])
        print('   Found data:', foundKeys)
        if 'Grain' in foundKeys:  # open new array if data exists
            self.grainID = -np.ones_like(self.phaseID)

        # read data
        data = np.loadtxt(fileName)
        print('   Reading file of size ', data.shape, '  this can take a bit...')
        if update:
            self.mask[:] = False
            print('Warning: this is too slow')
            """
      be intelligent where you seearch, check if old and new data monotonically increases
      then search in sections of equal y
      or subdivide into half, of half of half
      """
            xs, ys = self.xy()
            for i in range(data.shape[0]):
                x = data[i, foundKeys['x,'] - 1].astype(float)
                y = data[i, foundKeys['x,'] - 0].astype(float)
                # identify index: nice and much much slowes
                idx = np.argmax(np.logical_and(np.abs(xs-x) < self.stepSizeX/10.0,
                                               # very save error of 10th of STEPSIZE
                                               np.abs(ys-y) < self.stepSizeX/10.0))
                # update
                self.mask[idx] = True
                self.phi1[idx] = data[i,
                                      foundKeys['phi1,'] - 1].astype(np.float16)
                self.phi[idx] = data[i, foundKeys['phi1,'] -
                                     0].astype(np.float16)
                self.phi2[idx] = data[i,
                                      foundKeys['phi1,'] + 1].astype(np.float16)
                if 'IQ' in foundKeys:
                    self.iq[idx] = data[i, foundKeys['IQ'] -
                                        1].astype(np.float16)
                if 'CI' in foundKeys:
                    self.ci[idx] = data[i, foundKeys['CI'] -
                                        1].astype(np.float16)
                if 'Fit' in foundKeys:
                    self.fit[idx] = data[i, foundKeys['Fit'] -
                                         1].astype(np.float16)
                if 'Phase' in foundKeys:
                    self.phaseID[idx] = data[i,
                                             foundKeys['Phase'] - 1].astype(np.float16)
                if 'sem' in foundKeys:
                    self.semSignal[idx] = data[i,
                                               foundKeys['sem'] - 1].astype(np.float16)
                if 'Grain' in foundKeys:
                    self.grainID[idx] = data[i,
                                             foundKeys['Grain'] - 1].astype(np.float16)
                # stepSizeX, width, height etc do not change
        else:  # read new
            self.phi1 = data[:, foundKeys['phi1,'] - 1].astype(np.float16)
            self.phi = data[:, foundKeys['phi1,'] - 0].astype(np.float16)
            self.phi2 = data[:, foundKeys['phi1,'] + 1].astype(np.float16)
            x = data[:, foundKeys['x,'] - 1].astype(float)
            y = data[:, foundKeys['x,'] - 0].astype(float)
            if 'IQ' in foundKeys:
                self.iq = data[:, foundKeys['IQ'] - 1].astype(np.float16)
            if 'CI' in foundKeys:
                self.ci = data[:, foundKeys['CI'] - 1].astype(np.float16)
            if 'Fit' in foundKeys:
                self.fit = data[:, foundKeys['Fit'] - 1].astype(np.float16)
            if 'Phase' in foundKeys:
                self.phaseID = data[:,foundKeys['Phase'] - 1].astype(np.uint8)
                self.phaseID += not self.phaseID.any()  # single-phase EDAX files use 0
            if 'sem' in foundKeys:
                self.semSignal = data[:,
                                      foundKeys['sem'] - 1].astype(np.float16)
            self.mask   = np.ones_like(x, dtype=bool)
            self.width  = max(x)
            self.height = max(y)
            self.ratio  = self.width/self.height
            self._setGrid(x, y)
        fileHandle.close()
        print('Duration loadTXT: ', int(np.round(time.time()-startTime)), 'sec')
        return


    def loadOSC(self, fileName: str = '') -> None:
        """
        Load .osc file; filename saved in self. No need to use it.
        Copied from mtex and translated into python
        Warning: SEMsignal not parsed correctly

        Args:
           fileName: file to read [default: self.fileName]
        """
        print('TODO: Symmetry has to be read and used')
        if fileName:
            self.fileName = fileName
        print('Load .osc file: ', self.fileName)

        # OSC stores its numeric values as little-endian 32-bit values.  Using
        # NumPy's platform-sized ``float`` (normally float64) desynchronizes the
        # reader after the first step-size field.
        startBytes = bytes.fromhex('B9 0B EF FF 02 00 00 00')
        raw = Path(self.fileName).read_bytes()
        header = np.frombuffer(raw, dtype='<u4', count=8)
        n = int(header[6])  # number of data points
        startPosition = raw.find(startBytes)
        if startPosition < 0:
            raise ValueError('OSC data-block marker was not found.')

        dataOffset = startPosition + len(startBytes)
        dataSize = n * 10 * np.dtype('<f4').itemsize
        remaining = len(raw) - dataOffset
        if remaining == dataSize + 8:
            # Current format: x and y step sizes directly precede the records.
            pass
        elif remaining == dataSize + 12:
            # Older format: a uint32 record-size/count field precedes them.
            dataOffset += 4
        else:
            raise ValueError(
                f'Unexpected OSC data-block size: expected {dataSize + 8} or '
                f'{dataSize + 12} bytes after the marker, found {remaining}.')

        self.stepSizeX, self.stepSizeY = np.frombuffer(
            raw, dtype='<f4', count=2, offset=dataOffset).astype(float)
        data = np.frombuffer(raw, dtype='<f4', count=n*10,
                             offset=dataOffset + 8).reshape(n, 10)
        self.phi1 = data[:, 0].astype(np.float16)
        self.phi = data[:, 1].astype(np.float16)
        self.phi2 = data[:, 2].astype(np.float16)
        self.iq = data[:, 5].astype(np.float16)
        self.ci = data[:, 6].astype(np.float16)
        self.phaseID = data[:, 7].astype(np.uint8)
        self.phaseID += not self.phaseID.any()  # single-phase EDAX files use 0
        self.semSignal = data[:, 8].astype(np.float16)  # SEMSignal
        self.fit = data[:, 9].astype(np.float16)  # Fit
        self.width  = float(max(data[:, 3]))
        self.height = float(max(data[:, 4]))
        self.ratio = self.width/self.height
        self._setGrid(data[:, 3].astype(float), data[:, 4].astype(float))
        del data
        return


    def loadCRC(self, fileName: str = '') -> None:
        """
        Load .crc file; filename saved in self. No need to use it.
        Copied from mtex and translated into python

        Args:
           fileName: file to read [default: self.fileName]; the .cpr file of the same name holds the metadata
        """
        if fileName:
            self.fileName = fileName
        cprFileName = self.fileName[:-4]+'.cpr'
        print('Load .crc file: ', self.fileName, cprFileName)
        if not os.path.exists(cprFileName):
            print('CPR file does not exist')
        cprFile = open(cprFileName)
        cprData: dict[str, dict[str, Any]] = {}
        for line in cprFile:
            line = line.strip()
            if line[0] == '[':
                title = line[1:-1].lower()
                cprData[title] = {}
                continue
            key, value = line.split('=')[0], line.split('=')[1]
            try:
                cprData[title][key.lower()] = float(value)
            except ValueError:
                cprData[title][key.lower()] = value.lower()
        cprFile.close()
        # print "META DATA",cprData
        self.stepSizeX = np.double(cprData['job']['griddistx'])
        self.stepSizeY = np.double(cprData['job']['griddisty'])
        xcells = int(cprData['job']['xcells'])
        ycells = int(cprData['job']['ycells'])
        numDataPoints = xcells * ycells
        self.width = xcells * self.stepSizeX
        self.height = ycells * self.stepSizeY
        self.ratio = self.width/self.height
        self.sym = [Symmetry(OXFORD_LAUE_GROUPS.get(cprData[f'phase{i}']['lauegroup'], ''))
                    for i in range(1, int(cprData['phases']['count'])+1)]

        # verify that data in correct order
        allColumnNames = [
            'X',                  # 1    4 bytes
            'Y',                  # 2       "
            'phi1',               # 3       "
            'Phi',                # 4       "
            'phi2',               # 5       "
            'MAD',                # 6       "
            'BC',                 # 7    1 byte
            'BS',                 # 8       "
            'Unknown',            # 9       "
            'Bands',              # 10      "
            'Error',              # 11      "
            'ReliabilityIndex']    # 12      "
        allDataType = np.ones((12,), dtype=int)
        allDataType[:6] = 4
        allDataType[-1] = 4
        columnNames, columnType = ['Phase'], [1]
        for k in range(int(cprData['fields']['count'])):
            order = int(cprData['fields']['field'+str(k+1)])-1
            if order <= 12:
                columnNames.append(allColumnNames[order])
                columnType.append(allDataType[order])
            else:
                columnNames.append('Unknown'+str(order))
                columnType.append(4)
        expectedColumns = [
            'Phase', 'phi1', 'Phi', 'phi2', 'MAD', 'BC', 'BS', 'Bands',
            'Error', 'ReliabilityIndex',
        ]
        if columnNames == expectedColumns:
            print('  CRC-Data in correct order')
        else:
            print('  WARNING! CRC-Data not in correct order! WARNING')
            print(
                f'    should be {expectedColumns}'
            )
            print('    is       ', columnNames)
            print('    if data missing at end, no problem')
        # print columnType

        # coordinates
        xCoordinates = np.arange(xcells)*self.stepSizeX
        yCoordinates = np.arange(ycells)*self.stepSizeY
        x, y = np.meshgrid(xCoordinates, yCoordinates)
        self._setGrid(x.flatten(), y.flatten())

        # read data from crcFile: packed little-endian records
        recordType = [('phase', 'u1'), ('phi1', '<f4'), ('phi', '<f4'), ('phi2', '<f4'), ('ci', '<f4'),
                      ('bc', 'u1'), ('bs', 'u1'), ('bands', 'u1'), ('error', 'u1')]
        if 'ReliabilityIndex' in columnNames:
            recordType.append(('ri', '<f4'))
        data = np.fromfile(self.fileName, dtype=recordType, count=numDataPoints)
        self.phaseID = data['phase'].copy()
        self.bc, self.bs, self.bands, self.error = (data[i].copy() for i in ('bc', 'bs', 'bands', 'error'))
        self.phi1, self.phi, self.phi2, self.ci = (data[i].astype(float) for i in ('phi1', 'phi', 'phi2', 'ci'))
        self.ri = data['ri'].astype(float) if 'ri' in (data.dtype.names or ()) else np.zeros(numDataPoints)
        self.iq, self.semSignal, self.fit = (np.zeros(numDataPoints) for _ in range(3))
        if np.max(self.phaseID) > len(self.sym):
            print('ERRRO in reading CRC: symmetries do not match', len(self.sym), np.max(self.phaseID))
        return


    def loadVoid(self, rotation: str) -> None:
        """
        rotation angles in degree

        Args:
           rotation: Euler angles "phi1|Phi|phi2" in degrees, optionally followed by "|spread" (standard deviation
                     of the random scatter in radians) and "|numberPerAxis" [default: 6]; without "|": no rotation
        """
        numPerAxis, distrib = 6, 0.0
        if '|' in rotation:
            values = [float(i) for i in rotation.split('|')]
            if len(values) == 3:
                phi1, phi, phi2 = np.radians(values)
            elif len(values) == 4:
                phi1, phi, phi2 = np.radians(values[:3])
                distrib = values[-1]
            elif len(values) == 5:
                phi1, phi, phi2 = np.radians(values[:3])
                distrib, numPerAxis = values[-2], int(values[-1])
            else:
                print('ERROR')
                return
            print('   Euler angles:', np.round(phi1, 2), np.round(phi, 2), np.round(phi2, 2),
                  '| distribution:', distrib, '| numberPerAxis:', numPerAxis)
        else:
            phi1, phi, phi2 = 0, 0, 0
        if distrib < 0.001:
            distrib = 0.001
        self.sym.append(Symmetry('cubic'))
        self.stepSizeX = 1.
        numDataPoints = int(numPerAxis**2)
        coordinates = np.arange(numPerAxis)*self.stepSizeX
        x, y = np.meshgrid(coordinates, coordinates)
        self._setGrid(x.flatten(), y.flatten())
        self.phaseID = np.ones((numDataPoints), dtype=np.uint8)
        self.phi1 = np.zeros((numDataPoints), dtype=float)+phi1 + \
            np.random.normal(loc=0, scale=distrib, size=numDataPoints)
        self.phi = np.zeros((numDataPoints), dtype=float)+phi + \
            np.random.normal(loc=0, scale=distrib, size=numDataPoints)
        self.phi2 = np.zeros((numDataPoints), dtype=float)+phi2 + \
            np.random.normal(loc=0, scale=distrib, size=numDataPoints)
        self.ci = np.ones((numDataPoints), dtype=float)
        self.width = np.max(x)
        self.height = np.max(y)
        self.ratio = self.width/self.height
        return

    # @}


    ##
    # @name Mask and Path routines
    # Masked areas are plotted in black. Initially no point is part of the mask;
    # all points are false.
    # @{
    def maskCI(self, ci: float) -> None:
        """
        masked all points off, which have a CI less than: good points=False, bad points=True

        Args:
           ci: critical CI
        """
        self.mask = self.ci > ci
        return


    def maskReset(self) -> None:
        """
        reset mask
        """
        self.mask = self.ci > -1
        return


    def removePointsOutsideMask(self) -> None:
        """
        Set points outside the mask invalid so they are read as nonexistent by
        OIM after export.
        """
        self.ci[~self.mask] = -1.0
        self.fit[~self.mask] = 180.0
        return


    def setVMask(self, every: int = 1) -> None:
        """
        fast preview, does not influence results in any way: maps show only every k-th row and column
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
        """
        crop visible area

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
        """
        x- and y-coordinates of points, computed from the grid parameters

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


    def calcKAM(self, layers: int = 1) -> None:
        """
        calculate Kerner Average Misorientation in DEGREES (because user focused)

        Args:
           layers: number of neighboring layers used for KAM (more: slower)
        """
        startTime = time.time()
        neighbors = self.neighbors()
        assert neighbors is not None
        angles = np.full(neighbors.shape, np.nan)
        for phase, sym in enumerate(self.sym):
            points = np.flatnonzero(self.phaseID == phase)
            if not sym.lattice or not len(points):
                continue
            symQ = sym.symmetryQuats()
            for iNeighbor in range(neighbors.shape[1]):
                misQ = self.quaternions[points].inv() * self.quaternions[neighbors[points, iNeighbor]]
                angles[points, iNeighbor] = np.min([(misQ*q).magnitude() for q in symQ], axis=0)
        # -10 would index points at the end of the map
        angles[(neighbors < 0) | (self.phaseID[neighbors] != self.phaseID[:, None])] = np.nan
        angles[self.ci[neighbors] == -1.0] = np.nan
        self.kam = np.degrees(np.nanmean(angles, axis=1))
        self.kam[self.ci == -1.0] = np.nan
        print('Duration KAM evaluation: ', int(
            np.round(time.time()-startTime)), 'sec')
        return


    def neighbors(self, idx: int | None = None, layers: int = 1) -> np.ndarray | None:
        """
        identify neighboring indexes

        Args:
           idx: index to find [if None: calculate all]
           layers: number of neighboring layers

        Returns:
         array of neighbors; invalid points have a value=-10
        """
        if layers != 1:
            print('number of layers not implemented')
            return None
        i = np.arange(self.nPoints) if idx is None else np.atleast_1d(idx)
        period = self.nColsOdd + self.nColsEven
        pair, j = np.divmod(i, period)
        even = (j >= self.nColsOdd).astype(int)  # 2nd row of the pair of rows
        row, col = 2*pair + even, j - even*self.nColsOdd
        if self.grid == 'HexGrid':
            # first of the two touching columns in the rows above and below; depends on shift direction
            shift = int(self.xOffset > 0)
            lower = col + np.where(even, shift-1, -shift)
            nRow = row[:, None] + np.array([-1, -1, 0, 0, 1, 1])
            nCol = np.stack([lower, lower+1, col-1, col+1, lower, lower+1], axis=1)
        else:
            nRow = row[:, None] + np.array([-1, 0, 0, 1])
            nCol = np.stack([col, col-1, col+1, col], axis=1)
        rowLength = np.where(nRow % 2, self.nColsEven, self.nColsOdd)
        neighbors = (nRow//2)*period + (nRow % 2)*self.nColsOdd + nCol
        valid = (nRow >= 0) & (nRow < self.nRows) & (nCol >= 0) & (nCol < rowLength) & (neighbors < self.nPoints)
        neighbors[~valid] = -10
        if idx is not None:
            return neighbors[0]
        return neighbors

    # @}


    def writeANG(self, fileName: str) -> None:
        """
        write body of ang file

        Args:
           fileName: file name
        """
        startTime = time.time()
        fileOut = open(fileName, 'w')
        fileOut.write('# MaterialName void\n')
        fileOut.write('# Formula \n')
        # adopt for HCP (fcc and bcc the same)
        fileOut.write('# Symmetry 43\n')
        fileOut.write('# LatticeConstants 1.0 1.0 1.0 90.0 90.0 90.0\n')
        fileOut.write('# NumberFamilies 4\n')
        fileOut.write('# khlFamilies 1 1 1 1 0.0\n')  # adopt for HCP
        fileOut.write('# khlFamilies 2 0 0 1 0.0\n')
        fileOut.write('# khlFamilies 2 2 0 1 0.0\n')
        fileOut.write('# khlFamilies 3 1 1 1 0.0\n')
        fileOut.write(f'#\n# GRID: {self.grid}\n#\n')
        xs, ys = self.xy()
        phaseID = self.phaseID if self.phaseID.max() > 1 else np.zeros_like(self.phaseID)
        for i in range(self.nPoints):
            phi1, phi, phi2 = tuple(asBungeEulers(self.quaternions[i]))
            fileOut.write(
                f' {phi1:8.5f} {phi:8.5f} {phi2:8.5f} {xs[i]:12.5f}'
                f' {ys[i]:12.5f} {self.iq[i]:8.3f} {self.ci[i]:6.3f}'
                f' {phaseID[i]:2d} {self.semSignal[i]:6d} {self.fit[i]:7.3f}\n'
            )
        fileOut.close()
        print('Duration writeANG: ', int(
            np.round(time.time()-startTime)), 'sec')
        return


    @property
    def x(self) -> np.ndarray:
        """
        x-coordinates of all points

        Returns:
           x-coordinates in [um]
        """
        return self.xy()[0]


    @property
    def y(self) -> np.ndarray:
        """
        y-coordinates of all points

        Returns:
           y-coordinates in [um]
        """
        return self.xy()[1]


    def _setGrid(self, x: np.ndarray, y: np.ndarray) -> None:
        """
        Derive grid parameters from row-major coordinates; the coordinates themselves are not stored.

        Args:
           x: x-coordinates of all points
           y: y-coordinates of all points
        """
        rowStarts = np.flatnonzero(np.diff(x) < 0) + 1
        if len(rowStarts) == 0:
            raise ValueError('EBSD data must contain more than one scan row.')
        self.nPoints = len(x)
        self.nRows = len(rowStarts) + 1
        self.nColsOdd = int(rowStarts[0])
        self.nColsEven = int(rowStarts[1] if len(rowStarts) > 1 else len(x)) - self.nColsOdd
        self.x0, self.y0 = float(x[0]), float(y[0])
        self.stepSizeX = float(x[1] - x[0])
        self.stepSizeY = float(y[self.nColsOdd] - y[0])
        self.xOffset = float(x[self.nColsOdd] - x[0])
        self.grid = 'HexGrid' if abs(self.xOffset) > self.stepSizeX/10 else 'SqrGrid'
        xGrid, yGrid = self.xy()
        if not (np.allclose(xGrid, x, atol=self.stepSizeX/10) and np.allclose(yGrid, y, atol=self.stepSizeY/10)):
            raise ValueError('EBSD data is not a complete, row-major rectangular or hexagonal grid.')
        return


    def _image(self, values: np.ndarray, widthPixel: int | None = None,
               interpolationType: str = 'nearest') -> tuple[np.ndarray, tuple[float, float, float, float]]:
        """
        Sample per-point values onto the pixels of an image of the visible area

        Args:
           values: one value (or one row of values) per point
           widthPixel: horizontal size of the image [default: one pixel per step; half a step for hexagonal grids]
           interpolationType: "nearest" uses the grid directly; others interpolate with scipy's griddata

        Returns:
           image of shape (height, width, ...), extent for imshow in [um]
        """
        xVisible, yVisible = self.xy(np.flatnonzero(self.vMask))
        xLo, xHi = xVisible.min()-self.stepSizeX/2, xVisible.max()+self.stepSizeX/2
        yLo, yHi = yVisible.min()-self.stepSizeY/2, yVisible.max()+self.stepSizeY/2
        self.ratio = (xHi-xLo) / (yHi-yLo)
        if widthPixel is None:
            pitch = self.stepSizeX/2 if self.grid == 'HexGrid' else self.stepSizeX
            widthPixel = int(round((xHi-xLo)/pitch))
            heightPixel = int(round((yHi-yLo)/self.stepSizeY))
        else:
            heightPixel = int(round(widthPixel/self.ratio))
        xAxis = xLo + (np.arange(widthPixel)+0.5)*(xHi-xLo)/widthPixel
        yAxis = yLo + (np.arange(heightPixel)+0.5)*(yHi-yLo)/heightPixel
        extent = (xLo, xHi, yHi, yLo)
        if interpolationType != 'nearest':
            x, y = np.meshgrid(xAxis, yAxis)
            return griddata(np.vstack((xVisible, yVisible)).T, values[self.vMask], (x, y), interpolationType), extent
        if self.vMaskEvery > 1:  # preview: pixels sample only every k-th row and column
            xBlock, yBlock = self.stepSizeX*self.vMaskEvery, self.stepSizeY*self.vMaskEvery
            xAxis = xLo + (np.floor((xAxis-xLo)/xBlock)+0.5)*xBlock
            yAxis = yLo + (np.floor((yAxis-yLo)/yBlock)+0.5)*yBlock
        return values[self._nearestIndex(xAxis[None, :], yAxis[:, None])], extent


    def _nearestIndex(self, x: Any, y: Any) -> np.ndarray:
        """
        Index of the grid point closest to x, y: row from y, then column within that row from x

        Args:
           x: x-coordinate(s) in [um]; broadcast with y
           y: y-coordinate(s) in [um]

        Returns:
           point indices, of the broadcast shape of x and y
        """
        row = np.clip(np.rint((np.asarray(y)-self.y0)/self.stepSizeY).astype(int), 0, self.nRows-1)
        even = row % 2
        col = np.rint((np.asarray(x) - self.x0 - even*self.xOffset)/self.stepSizeX).astype(int)
        col = np.clip(col, 0, np.where(even, self.nColsEven, self.nColsOdd)-1)
        idx = (row//2)*(self.nColsOdd+self.nColsEven) + even*self.nColsOdd + col
        return np.minimum(idx, self.nPoints-1)
