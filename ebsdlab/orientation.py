##
# @file
# @brief Orientation class: combination of material symmetry and specific rotation.
# Copyright:
#   Original version was part of DAMASK <damask.mpie.de> (Martin Diehl, Philip Eisenlohr, Franz Roters)
#
import numpy as np
from scipy.spatial.transform import Rotation
from ._rotation import asBungeEulers
from .symmetry import Symmetry


class Orientation:
    """Orientation class: combination of material symmetry and specific rotation
    """
    __slots__ = ['quaternion', 'symmetry', 'plot2D', 'eps']

    # @name CONVENTIONAL ROUTINES
    # @{

    def __init__(self,
                 quaternion=None,
                 matrix=None,
                 eulers=None,
                 # put any integer to have a fixed seed or True for real random
                 random=False,
                 symmetry='',
                 ):
        # produce random orientation
        if random:
            if isinstance(random, bool):
                self.quaternion = Rotation.random()
            else:
                self.quaternion = Rotation.random(random_state=random)
        # based on given Euler angles
        elif isinstance(eulers, np.ndarray) and eulers.shape == (3,):
            self.quaternion = Rotation.from_euler('ZXZ', eulers)
        # based on given rotation matrix
        elif isinstance(matrix, np.ndarray):
            self.quaternion = Rotation.from_matrix(matrix)
        # based on a SciPy rotation
        elif isinstance(quaternion, Rotation):
            self.quaternion = Rotation.from_quat(quaternion.as_quat())
        else:
            self.quaternion = Rotation.identity()
        self.symmetry = Symmetry(symmetry)
        self.plot2D = 'down-right'
        self.eps = 1e-6
        return

    def __copy__(self):
        return self.__class__(quaternion=self.quaternion, symmetry=self.symmetry.lattice)

    def __repr__(self):
        matrix = '\n'.join('\t'.join(map(str, self.asMatrix()[i, :])) for i in range(3))
        eulers = '\t'.join(map(str, self.asEulers('bunge', degrees=True)))
        return (f'Symmetry: {self.symmetry}\n'
                f'Quaternion: {self.quaternion}\n'
                f'Matrix:\n{matrix}\n'
                f'Bunge Eulers / deg: {eulers}')

    def asEulers(self,
                 notation='bunge',
                 degrees=False,
                 standardRange=False):
        """Return this orientation's Euler angles in the requested convention."""
        if notation.lower() not in ('bunge', 'zxz'):
            raise ValueError("Only the Bunge/intrinsic ZXZ convention is supported")
        return asBungeEulers(self.quaternion, degrees, standardRange)
    eulers = property(asEulers)

    def asMatrix(self):
        """Return this orientation as a 3-by-3 rotation matrix."""
        return self.quaternion.as_matrix()
    matrix = property(asMatrix)

    def inFZ(self):
        """Check whether given Rodrigues vector falls into fundamental zone of own symmetry.
        """
        return self.symmetry.inFZ(self.quaternion)
    infz = property(inFZ)

    def equivalentQuaternions(self, who=None):
        """Return symmetry-equivalent quaternions, optionally selected by index."""
        return self.symmetry.equivalentQuaternions(self.quaternion, who)

    def equivalentOrientations(self, who=None):
        """Return symmetry-equivalent orientations, optionally selected by index."""
        return [
            Orientation(quaternion=q, symmetry=self.symmetry.lattice)
            for q in self.equivalentQuaternions(who)
        ]

    def reduced(self):
        '''
        Transform orientation to fall into fundamental zone according to symmetry
        '''
        for me in self.symmetry.equivalentQuaternions(self.quaternion):
            if self.symmetry.inFZ(me):
                break
        return Orientation(quaternion=me, symmetry=self.symmetry.lattice)

    # @}
    ##
    # @name MATERIAL SPECIFIC ROUTINES
    # @{

    def disorientation(self, other, sst=True):
        """Disorientation between myself and given other orientation.

        Rotation axis falls into SST if SST == True.
          (Currently requires same symmetry for both orientations.
          Look into A. Heinz and P. Neumann 1991 for cases with differing sym.)

         Args:
          other: other orientation
          SST: True (rotation axis falls into SST); False

        Returns:
          disorientation quaternion; indices of equivalent orientations; and
          whether the result was conjugated
        """
        if self.symmetry != other.symmetry:
            raise TypeError(
                'disorientation between different symmetry classes not supported yet.')
        misQ = self.quaternion.inv()*other.quaternion
        mySymQs = self.symmetry.symmetryQuats() if sst else self.symmetry.symmetryQuats()[
            :1]       # take all or only first sym operation
        otherSymQs = other.symmetry.symmetryQuats()
        for i, sA in enumerate(mySymQs):  # if not in SST: only one sA
            for j, sB in enumerate(otherSymQs):  # changes always
                candidate = sA.inv()*misQ*sB
                for k, theQ in enumerate((candidate.inv(), candidate)):
                    breaker = self.symmetry.inFZ(theQ) and (
                        not sst or other.symmetry.inDisorientationSST(theQ))
                    if breaker:
                        break
                if breaker:
                    break
            if breaker:
                break
        return (Orientation(quaternion=theQ, symmetry=self.symmetry.lattice),
                # disorientation, own sym, other sym, self-->other: True, self<--other: False
                i, j, k == 1)

    def inversePole(self, axis, proper=False, sst=True):
        """Rotate an axis into the standard stereographic triangle using symmetry.

        Args:
          axis: vector in crystal orientation, e.g. [100]
          proper: consider only vectors with z >= 0 using two neighboring SSTs;
              this permits more positive results without changing the RGB value
          SST: iterate through all equivalent and find the one in the SST

        Returns:
          vector of axis
        """
        if sst:  # Pole requested to be within SST.
            # test all symmetric equivalent quaternions
            for i, q in enumerate(self.symmetry.equivalentQuaternions(self.quaternion)):
                # align crystal direction to axis
                pole = q.inv().apply(axis)
                if self.symmetry.inSST(pole, proper):
                    break                                                # found SST version
        else:
            # align crystal direction to axis
            pole = self.quaternion.inv().apply(axis)
        return (pole, i if sst else 0)

    def ipfColor(self, axis, proper=False):
        """color of inverse pole figure for given axis

        Args:
           axis: axis of pole figure (ND=001)
           proper: consider vectors with z >= 0 using two neighboring SSTs;
               this permits more positive results without changing the RGB value

        Returns:
           vector of color (rgb)
        """
        color = np.zeros(3, 'd')
        for q in self.symmetry.equivalentQuaternions(self.quaternion):
            # align crystal direction to axis
            pole = q.inv().apply(axis)
            inSST, color = self.symmetry.inSST(pole, color=True, proper=proper)
            if inSST:
                break
        return color

    @classmethod
    def average(cls,
                orientations,
                multiplicity=None):
        """Return the average orientation

        ref: F. Landis Markley, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman,
          Averaging Quaternions,
          Journal of Guidance, Control, and Dynamics, Vol. 30, No. 4 (2007), pp. 1193-1197.
          doi: 10.2514/1.28949

        Usage:
          * a = Orientation(eulers=np.radians([10, 10, 0]), symmetry='hexagonal')
          * b = Orientation(eulers=np.radians([20, 0, 0]),  symmetry='hexagonal')
          * avg = Orientation.average([a,b])

        Args:
          cls: class method (void)
          orientations: list of orientations
          multiplicity: --

        Returns:
          average orientation (not rotation) in radians
        """
        if not all(isinstance(item, Orientation) for item in orientations):
            raise TypeError('Only instances of Orientation can be averaged.')
        count = len(orientations)
        if multiplicity is None or len(multiplicity) == 0:
            multiplicity = np.ones(count, dtype='i')
        # take first as reference
        reference = orientations[0]
        closestRotations = []
        for o in orientations:
            closest = o.equivalentOrientations(reference.disorientation(o, sst=False)[2])[
                0]             # select sym orientation with lowest misorientation
            closestRotations.append(closest.quaternion)
        mean = Rotation.concatenate(closestRotations).mean(weights=multiplicity)
        return Orientation(quaternion=mean,
                           symmetry=reference.symmetry.lattice)

    # @}
    ##
    # @name PLOTTING, PRINTING
    # @{
    def project(self, x, y, z):
        """

        down-right: y, -x
        up-left   : -y, x
        right-up   : x,y
        left-down: -x,-y
        3D        : x,y,z
        """
        x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)
        projections = {'down-right': (y, -x), 'up-left': (-y, x), 'right-up': (x, y),
                       'left-down': (-x, -y), '3D': (x, y, z)}
        if self.plot2D not in projections:
            print('Error: plot2D not well defined: plotLine')
            return None
        return projections[self.plot2D]

    def plotLine(self, ax, start, delta, color='k', lw=1, ls='solid', markerSize=None):
        """
        Plot one line using given projection

        Args:
           ax: axis to plot into
           start: start coordinate
           delta: delta cooradinate (end-start)
           color: color
           lw: line width
           ls: line style "solid",'dashed'
           markerSize: size of marker (only used for non-lines: delta>0)
        """
        if np.linalg.norm(delta) < self.eps:
            marker = 'o'
            if markerSize is None:
                markerSize = 7
        else:
            marker = None
            markerSize = 0
        if self.plot2D == '3D':
            ax.plot([start[0]]+[start[0]+delta[0]],
                    [start[1]]+[start[1]+delta[1]],
                    [start[2]]+[start[2]+delta[2]],
                    color=color, lw=lw, marker=marker, ls=ls, markersize=markerSize)
        else:
            x, y = self.project([start[0]]+[start[0]+delta[0]],
                                [start[1]]+[start[1]+delta[1]],
                                [start[2]]+[start[2]+delta[2]])
            ax.plot(x, y,  color=color, lw=lw, marker=marker,
                    ls=ls, markersize=markerSize)
        return

    def plotUnit(self, ax, xlabel, ylabel, zlabel,  x=0, y=0, z=0,   s=1):  # unit axis
        """
        Coordinate systems: see plotLine

        Args:
           ax: axis to used for plotting
           xlabel: x-label
           ylabel: y-label
           zlabel: z-label
           x: x-coordinate of origin
           y: y-coordinate of origin
           z: z-coordinate of origin
           s: scale
        """
        self.plotLine(ax, [x, y, z], [s, 0, 0], 'k', lw=3)
        self.plotLine(ax, [x, y, z], [0, s, 0], 'k', lw=3)
        self.plotLine(ax, [x, y, z], [0, 0, s], 'k', lw=3)
        if self.plot2D == '3D':
            ax.text(x+s,   y+0.1, z+0.1, xlabel)
            ax.text(x+0.1, y+s,   z+0.1, ylabel)
            ax.text(x+0.1, y+0.1, z+s, zlabel)
        else:
            ax.text(*(self.project(x+s,   y, z+0.1)+(xlabel,)))
            ax.text(*(self.project(x+0.1, y+s, z+0.1) +
                    (ylabel, {'ha': 'right'})))
            ax.text(*(self.project(x+0.1, y, z+s)+(zlabel,)))
        return

    def plot(self, poles=None, unitCell=True, cos=True, annotate=False, plot2D='', scale=2, fileName=''):
        """Plot rotated unit-cell in 3D, and possibly the pole-figure and specific poles

        Projection onto 2D: cooradinate systems are given as xDirection-yDirection (z follows)
        - down-right: [default in text books] RD = x = down; TD = y = right; ND = z = outOfPlane
        - up-left: [default in OIM] RD = x = up; TD = y = left; ND = z = outOfPlane

        Args:
           poles: if given (e.g. [1,0,0]), plot pole-figure and the corresponding poles
           unitCell: plot unit cell
           cos: plot coordinate system
           annotate: annotate poles in pole figure (requires poles given)
           plot2D: do a normal projection onto 2D plane: [down-right, up-left, right-up, left-down, 3D]; '' keeps the current setting
           scale: scale of pole-figure dome over crystal
           fileName: fileName for image output (if given, image not shown)
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        if plot2D:
            self.plot2D = plot2D
        if self.plot2D == '3D':
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            # ax.view_init(90,0)
        else:
            fig, ax = plt.subplots()

        if unitCell:
            for line in self.symmetry.unitCell():
                start, end = np.array(line[:3], dtype=float), np.array(
                    line[3:], dtype=float)
                start = self.quaternion.apply(start)
                end = self.quaternion.apply(end)
                if start[2] < 0 and end[2] < 0:
                    self.plotLine(ax, start, end-start, color='b', lw=0.2)
                elif start[2] > 0 and end[2] > 0:
                    self.plotLine(ax, start, end-start, color='b', lw=2)
                else:
                    delta = end-start
                    k = -start[2]/delta[2]
                    mid = start+k*delta
                    if start[2] > 0:
                        self.plotLine(ax, start, mid-start, color='b', lw=2)
                        self.plotLine(ax, mid,   end-mid, color='b', lw=0.2)
                    else:
                        self.plotLine(ax, start, mid-start, color='b', lw=0.2)
                        self.plotLine(ax, mid,   end-mid, color='b', lw=2)

        if cos:
            self.plotUnit(ax, 'RD [100]', 'TD [010]', 'ND [001]')

        if poles is not None:
            # plot sphere
            if self.plot2D == '3D':
                u = np.linspace(0, 2*np.pi, 50)
                v = np.linspace(0, np.pi/2, 50)
                x = scale*np.outer(np.cos(u), np.sin(v))
                y = scale*np.outer(np.sin(u), np.sin(v))
                z = scale*np.outer(np.ones_like(u), np.cos(v))
                ax.plot_surface(x, y, z, color='gray', alpha=0.7,
                                cstride=10, rstride=10, lw=0)
            else:
                ax.plot(scale*np.cos(np.linspace(0., 2.*np.pi, 100)),
                        scale*np.sin(np.linspace(0., 2.*np.pi, 100)), 'k--')
            # plot poles
            oHelp = Orientation(eulers=np.array(
                [0., 0., 0.]), symmetry=self.symmetry.lattice)
            poles = np.array(poles, dtype=float)
            poles /= np.linalg.norm(poles)
            for _, q in enumerate(oHelp.symmetry.equivalentQuaternions(oHelp.quaternion)):
                conjAxis = q.apply(poles)  # e.g. [100]
                direction = self.quaternion.apply(conjAxis)
                if direction[2] < -self.eps:
                    continue  # prevent rounding errors
                fromBase = direction+np.array([0, 0, 1])
                # self.plotLine(ax, [0,0,0], direction*scale, color='c', lw=1) #in plane lines: not needed
                if self.plot2D == '3D':
                    # lines from bottom base to points
                    self.plotLine(
                        ax, -scale*np.array([0, 0, 1]), (fromBase)*scale, 'c')
                xy = fromBase/fromBase[2]*scale
                xy[2] = 0.0
                self.plotLine(ax, xy, [0., 0., 0.], color='c')  # plot point
                if annotate:
                    xCoordinate, yCoordinate = self.project(xy[0], xy[1], xy[2])
                    label = str(np.array(conjAxis, dtype=int))[1:-1]
                    label = label.replace(' ', '')
                    ax.text(xCoordinate+0.05, yCoordinate+0.05, label)

        # finalize plot
        if self.plot2D == '3D':
            ax.axis('equal')
        else:
            ax.set_aspect('equal', adjustable='box')
        ax.axis('off')
        ax.set_xlim([-scale*1.1, scale*1.1])
        ax.set_ylim([-scale*1.1, scale*1.1])
        if self.plot2D == '3D':
            ax.set_zlabel('')
            ax.set_zticks([])
        if fileName:
            plt.savefig(fileName, dpi=150, bbox_inches='tight')
        else:
            plt.show()
        return

    def toScreen(self, equivalent=True):
        """
        print Euler angles and HKL /UVW

        Args:
          equivalent: print also equivalent orientations
        """
        print('Euler angles:', np.round(asBungeEulers(self.quaternion, degrees=True), 1))
        rotM = self.quaternion.as_matrix()
        print('HKL', np.array(rotM[2, :]/np.min(rotM[2, :]), dtype=int))
        print('UVW', -np.array(rotM[0, :]/np.min(rotM[0, :]), dtype=int))
        if equivalent:
            print('Equivalent orientations - Euler angles:')
            for q in self.symmetry.equivalentQuaternions(self.quaternion):
                angles = asBungeEulers(q, degrees=True)
                angles[angles < 0] += 360.
                print(f'   [{angles[0]:5.1f}  {angles[1]:5.1f}  {angles[2]:5.1f}]')
        return


    # @}
