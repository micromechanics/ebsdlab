"""Crystal symmetry operations and visualization helpers."""

import math
from typing import Any
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation
from ._rotation import asRodrigues

# Directions on an SST edge might give components of -1e-15 instead of exactly 0; count them as inside
SST_TOLERANCE = 1e-12

# Standard stereographic triangle (SST): smallest region representing
# symmetry-equivalent crystal directions in an inverse pole figure (IPF).
# The improper basis uses inversion/reflection equivalence and folds opposite directions into one triangle.
# The proper basis describes the adjoining triangle needed when equivalence is
# restricted to handedness-preserving rotations.
# Note: SciPy's rotation groups themselves contain proper rotations only.
GROUPS: dict[str, dict[str, Any]] = {
    'cubic': {
        'rotation_group': 'O',
        'sst_bases': {
            'improper': np.array([
                [-1.0, 0.0, 1.0],
                [np.sqrt(2.0), -np.sqrt(2.0), 0.0],
                [0.0, np.sqrt(3.0), 0.0],
            ]),
            'proper': np.array([
                [0.0, -1.0, 1.0],
                [-np.sqrt(2.0), np.sqrt(2.0), 0.0],
                [np.sqrt(3.0), 0.0, 0.0],
            ]),
        },
        'cell': {
            'geometry': 'orthogonal',
            'default_ratio': (1.0, 1.0, 1.0),
            'equal_to_a': ('b', 'c'),
            'default_angles': (90.0, 90.0, 90.0),
            'fixed_angles': (90.0, 90.0, 90.0),
        },
    },
    'hexagonal': {
        'rotation_group': 'D6',
        'sst_bases': {
            'improper': np.array([
                [0.0, 0.0, 1.0],
                [1.0, -np.sqrt(3.0), 0.0],
                [0.0, 2.0, 0.0],
            ]),
            'proper': np.array([
                [0.0, 0.0, 1.0],
                [-1.0, np.sqrt(3.0), 0.0],
                [np.sqrt(3.0), -1.0, 0.0],
            ]),
        },
        'cell': {
            'geometry': 'hexagonal',
            'default_ratio': (1.0, 1.0, 1.5),
            'equal_to_a': ('b',),
            'default_angles': (90.0, 90.0, 120.0),
            'fixed_angles': (90.0, 90.0, 120.0),
        },
    },
    'tetragonal': {
        'rotation_group': 'D4',
        'sst_bases': {
            'improper': np.array([
                [0.0, 0.0, 1.0],
                [1.0, -1.0, 0.0],
                [0.0, np.sqrt(2.0), 0.0],
            ]),
            'proper': np.array([
                [0.0, 0.0, 1.0],
                [-1.0, 1.0, 0.0],
                [np.sqrt(2.0), 0.0, 0.0],
            ]),
        },
        'cell': {
            'geometry': 'orthogonal',
            'default_ratio': (1.0, 1.0, 1.5),
            'equal_to_a': ('b',),
            'default_angles': (90.0, 90.0, 90.0),
            'fixed_angles': (90.0, 90.0, 90.0),
        },
    },
    'orthorhombic': {
        'rotation_group': 'D2',
        'sst_bases': {
            'improper': np.array([
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]),
            'proper': np.array([
                [0.0, 0.0, 1.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]),
        },
        'cell': {
            'geometry': 'orthogonal',
            'default_ratio': (1.0, 1.25, 1.5),
            'equal_to_a': (),
            'default_angles': (90.0, 90.0, 90.0),
            'fixed_angles': (90.0, 90.0, 90.0),
        },
    },
    'monoclinic': {
        'rotation_group': 'C2',
        'rotation_axis': 'Y',
        'sst_bases': {
            # Unique axis b is represented by the y axis.  The two regions
            # differ only by the sign of y and together form the proper FZ.
            'improper': np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                  [0.0, 1.0, 0.0]]),
            'proper': np.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0],
                                [0.0, -1.0, 0.0]]),
        },
        'cell': {
            'geometry': 'monoclinic',
            'default_ratio': (1.0, 1.25, 1.5),
            'equal_to_a': (),
            'default_angles': (90.0, 105.0, 90.0),
            'fixed_angles': (90.0, None, 90.0),
        },
    },
    'triclinic': {
        'rotation_group': 'C1',
        'sst_bases': {
            # With only inversion equivalence, an SST is one hemisphere.
            'improper': np.array([[0.0, 0.0, 1.0]] * 3),
            'proper': np.array([[0.0, 0.0, -1.0]] * 3),
        },
        'cell': {
            'geometry': 'triclinic',
            'default_ratio': (1.0, 1.25, 1.5),
            'equal_to_a': (),
            'default_angles': (75.0, 100.0, 110.0),
            'fixed_angles': (None, None, None),
        },
    },
    'trigonal': {
        'rotation_group': 'D3',
        'sst_bases': {
            # A 60 degree azimuthal wedge in the upper/lower hemisphere, between the mirror planes of D3 (2-fold
            # axes along a = x): from [10-10] at 30° (blue) to [01-10] at 90° (green), as Oxford's key 001-210-120
            'improper': np.array([[0.0, 0.0, 1.0], [-1.0, np.sqrt(3.0), 0.0], [2.0, 0.0, 0.0]]),
            'proper': np.array([[0.0, 0.0, -1.0], [-1.0, np.sqrt(3.0), 0.0], [2.0, 0.0, 0.0]]),
        },
        'cell': {
            'geometry': 'rhombohedral',
            'default_ratio': (1.0, 1.0, 1.0),
            'equal_to_a': ('b', 'c'),
            'default_angles': (75.0, 75.0, 75.0),
            'fixed_angles': (None, None, None),
        },
    },
}

# Rhombohedral is the lattice name customarily used for the trigonal crystal system.
LATTICE_ALIASES = {'rhombohedral': 'trigonal'}


class Symmetry:
    """Material symmetry identified by a human-readable lattice name."""

    def __init__(self, symmetry: str = '') -> None:
        """Create a symmetry from a lattice name such as 'cubic'; '' means not identified.

        Args:
           symmetry: lattice name or alias such as 'rhombohedral'
        """
        self.lattice = ''
        if symmetry:
            lattice = LATTICE_ALIASES.get(symmetry.lower(), symmetry.lower())
            if lattice not in GROUPS:
                raise ValueError(f'Unknown symmetry {symmetry!r}. Supported: {', '.join(GROUPS)}')
            self.lattice = lattice


    def standardTriangle(self, fileName: str = '', show: bool = True, stepSize: float = 0.1) -> Any:
        """Plot the colored cubic standard stereographic triangle.

        Args:
           fileName: if given, save to file
           show: show the figure
           stepSize: spacing of the colored points

        Returns:
           matplotlib figure; None if the lattice is not cubic
        """
        if self.lattice != 'cubic':
            print('ERROR: only implemented for cubic lattice')
            return None
        # border: [100]-[110] line and the arc from [110] to [111] (0 to 15 degrees)
        angles = np.radians(np.arange(16))
        arcX, arcY = np.cos(angles)*np.sqrt(2.0) - 1, np.sin(angles)*np.sqrt(2.0)
        border = np.vstack(([0.0, 0.0], np.column_stack((arcX, arcY)), [0.0, 0.0]))
        # colored background: grid points between the [100]-[111] line (x >= y) and the arc
        y, x = (a.ravel() for a in np.meshgrid(np.arange(0, 0.366025403784, stepSize),
                                               np.arange(0, 0.41421, stepSize), indexing='ij'))
        inside = (x >= y) & (x < np.interp(y, arcY, arcX))
        xy = np.stack((x[inside], y[inside]))
        hkl = self.xyToHKL(xy)
        assert hkl is not None
        _, rgb = self.inSST(hkl, color=True, proper=False)
        # plotting: background, border, labels
        fig, ax = plt.subplots()
        ax.scatter(xy[0], xy[1], c=rgb.T, s=15000. *stepSize, linewidths=0)
        ax.plot(border[:, 0], border[:, 1], '-k', linewidth=3)
        ax.text(0, 0, '[100]', horizontalalignment='right', fontsize=18)
        ax.text(0.42, 0, '[110]', fontsize=18)
        ax.text(0.37, 0.37, '[111]', fontsize=18)
        ax.axis('off')
        ax.axis('equal')
        showOrSave(fileName, show)
        return fig


    def symmetryQuats(self, who: Any = None) -> Rotation:
        """Return the proper symmetry rotations, with identity at index zero.

        Args:
           who: indices of the operations to return [default: all]

        Returns:
           symmetry rotations
        """
        # Generate all symmetry-equivalent rotation operations
        config = GROUPS.get(self.lattice)
        operations = (Rotation.create_group(config['rotation_group'], axis=config.get('rotation_axis', 'Z'))
                      if config is not None else Rotation.identity(1))
        # SciPy's cubic group does not start with identity, but callers rely on
        #    operations[0] being the rotation that leaves vectors unchanged.
        identity = int(np.argmin(operations.magnitude()))
        if identity != 0:
            order = np.concatenate(([identity], np.delete(np.arange(len(operations)), identity)))
            operations = operations[order]
        # Select a subset of those equivalent rotations
        if who is not None and np.size(who) > 0:
            operations = operations[np.atleast_1d(who)]
        return operations


    def unitCell(self, *, a: float = 1.0, b: float | None = None, c: float | None = None,
                 alpha: float | None = None, beta: float | None = None, gamma: float | None = None
                 ) -> np.ndarray | list[list[None]]:
        """Return centered unit-cell edges as ``[x1, y1, z1, x2, y2, z2]``.

        Args:
           a: lattice constant a
           b: lattice constant b [default: from the lattice's default ratio]
           c: lattice constant c [default: from the lattice's default ratio]
           alpha: lattice angle in degrees [default: lattice default]
           beta: lattice angle in degrees [default: lattice default]
           gamma: lattice angle in degrees [default: lattice default]

        Returns:
           edges, one row per edge; [[None]] for an unknown lattice
        """
        if self.lattice not in GROUPS:
            print('Unit cell not implemented')
            return [[None]]
        cellConfig = GROUPS[self.lattice]['cell']
        _, defaultB, defaultC = cellConfig['default_ratio']
        b = a * defaultB if b is None else b
        c = a * defaultC if c is None else c
        dimensions = np.asarray((a, b, c), dtype=float)
        if not np.all(np.isfinite(dimensions)) or np.any(dimensions <= 0.0):
            raise ValueError('lattice constants must be finite and positive')
        axisValues = {'b': b, 'c': c}
        for axis in cellConfig['equal_to_a']:
            if not np.isclose(a, axisValues[axis]):
                raise ValueError(f'{self.lattice} requires a = {axis}')
        defaultAngles = cellConfig.get('default_angles', (90.0, 90.0, 90.0))
        angles = np.asarray(tuple(default if value is None else value
                                  for value, default in zip((alpha, beta, gamma), defaultAngles)), dtype=float)
        if not np.all(np.isfinite(angles)) or np.any((angles <= 0.0) | (angles >= 180.0)):
            raise ValueError('lattice angles must be finite and between 0 and 180 degrees')
        for angle, required, name in zip(angles, cellConfig.get('fixed_angles', (None,) * 3),
                                         ('alpha', 'beta', 'gamma')):
            if required is not None and not np.isclose(angle, required):
                raise ValueError(f'{self.lattice} requires {name} = {required}')
        if cellConfig['geometry'] == 'rhombohedral' and not np.allclose(angles, angles[0]):
            raise ValueError('trigonal requires alpha = beta = gamma')
        if cellConfig['geometry'] == 'hexagonal':
            return _hexagonalCell(a, c)
        return _parallelepipedCell(a, b, c, *angles)


    def unitCellSegments(self, rotation: Rotation, scale: float = 1.) -> list[tuple[np.ndarray, np.ndarray, float]]:
        """Return the rotated unit-cell edges, split where they cross z=0.

        Args:
           rotation: orientation of the unit cell
           scale: scale of the unit cell

        Returns:
           segments as (start, end, line width): 2 above z=0, 0.2 below
        """
        segments = []
        for line in self.unitCell():
            start = rotation.apply(np.array(line[:3], dtype=float)*scale)
            end   = rotation.apply(np.array(line[3:], dtype=float)*scale)
            if start[2] < 0 and end[2] < 0:
                segments.append((start, end, 0.2))
            elif start[2] > 0 and end[2] > 0:
                segments.append((start, end, 2))
            else:
                delta = end-start
                mid = start+(-start[2]/delta[2])*delta
                if start[2] > 0:
                    segments += [(start, mid, 2), (mid, end, 0.2)]
                else:
                    segments += [(start, mid, 0.2), (mid, end, 2)]
        return segments


    def equivalentQuaternions(self, quaternion: Rotation, who: Any = None) -> list[Rotation]:
        """Return rotations equivalent to ``quaternion`` under this symmetry.

        Args:
           quaternion: rotation
           who: indices of the symmetry operations [default: all]

        Returns:
           equivalent rotations
        """
        return [quaternion*q for q in self.symmetryQuats(who)]


    def inFZ(self, rotation: Rotation) -> bool:
        """Return whether a rotation lies in the fundamental zone (FZ)

        Args:
           rotation: rotation to test

        Returns:
           True if inside the fundamental zone
        """
        # the FZ of a proper rotation group is the Voronoi region of identity: no equivalent rotation is smaller
        magnitudes = (rotation * self.symmetryQuats()).magnitude()
        return bool(magnitudes[0] <= np.min(magnitudes) + 1e-12)


    def inDisorientationSST(self, rotation: Rotation) -> bool:
        """Return whether a misorientation lies in the standard triangle.

        The criteria follow Heinz and Neumann, Acta Cryst. A47 (1991), 780-789.

        Args:
           rotation: misorientation to test

        Returns:
           True if inside the standard triangle
        """
        rodrigues = asRodrigues(rotation)
        if self.lattice == 'cubic':
            return bool(
                rodrigues[0] >= rodrigues[1]
                and rodrigues[1] >= rodrigues[2]
                and rodrigues[2] >= 0.0
            )
        if self.lattice == 'hexagonal':
            return bool(
                rodrigues[0] >= math.sqrt(3.0)*rodrigues[1]
                and rodrigues[1] >= 0.0 and rodrigues[2] >= 0.0
            )
        if self.lattice == 'tetragonal':
            return bool(
                rodrigues[0] >= rodrigues[1]
                and rodrigues[1] >= 0.0
                and rodrigues[2] >= 0.0
            )
        if self.lattice == 'orthorhombic':
            return bool(
                rodrigues[0] >= 0.0
                and rodrigues[1] >= 0.0
                and rodrigues[2] >= 0.0
            )
        return True


    def inSST(self, vector: Any, proper: bool = False, color: bool = False) -> Any:
        """Test vectors against the standard triangle and optionally color them.

        ``proper`` also considers the neighboring proper triangle. With
        ``color=True``, return the membership flags and IPF colors as RGB values.

        Args:
           vector: one vector (3,) or several vectors (3, n)
           proper: also consider the neighboring proper triangle
           color: also return IPF colors

        Returns:
           membership flag(s); with color: (flags, RGB colors of shape (3,) or (3, n))
        """
        config = GROUPS.get(self.lattice)
        basis = config['sst_bases'] if config is not None else None
        # theComponents = color components; inSST = membership flags
        inSST: Any
        if basis is None:
            theComponents = -np.ones(3, 'd')
            inSST = False
        else:
            v = np.array(vector, dtype=float)
            # check both improper ...
            if proper:
                theComponents = np.dot(basis['improper'], v)
                inSST = np.all(theComponents >= -SST_TOLERANCE, axis=0)
                # ... and proper SST
                if not np.any(inSST):
                    theComponents = np.dot(basis['proper'], v)
            else:
                theComponents = np.dot(basis['improper'], v)
                inSST = np.all(theComponents >= -SST_TOLERANCE, axis=0)
                oppositeComponents = np.dot(basis['improper'], -v)
                oppositeSST = np.all(oppositeComponents >= -SST_TOLERANCE, axis=0)
                if np.ndim(inSST) == 0:
                    if not inSST and oppositeSST:
                        theComponents = oppositeComponents
                else:
                    useOpposite = ~inSST & oppositeSST
                    theComponents[:, useOpposite] = oppositeComponents[:, useOpposite]
            inSST = np.all(theComponents >= -SST_TOLERANCE, axis=0)
        # have to return color array
        if color:
            if np.any(inSST):
                theComponentsNorm = theComponents / np.linalg.norm(theComponents, axis=0)
                # smoothen color ramps
                rgb = np.power(np.abs(theComponentsNorm), 0.5)
                if rgb.ndim > 1:
                    rgb[:, ~inSST] = 0.0
                if rgb.ndim > 1:
                    # normalize to (HS)V = 1
                    rgb[:, inSST] /= np.max(rgb[:, inSST], axis=0)
                else:
                    rgb /= np.max(rgb)
            else:
                rgb = np.zeros(3, 'd')
            return (inSST, rgb)
        return inSST


    def foldIntoSST(self, vectors: np.ndarray) -> np.ndarray:
        """Map directions onto their symmetry-equivalent direction in the standard triangle.

        Same result as trying all symmetryQuats() with inSST(proper=False), which includes inversion, but without
        the loop: the Laue group is applied directly by sign changes, sorting and an azimuthal fold.

        Args:
           vectors: directions (3, n)

        Returns:
           folded directions (3, n)
        """
        x, y, z = np.array(vectors, dtype=float)
        if self.lattice == 'cubic':  # 0 <= y <= x <= z
            y, x, z = np.sort(np.abs([x, y, z]), axis=0)
        elif self.lattice == 'tetragonal':  # 0 <= y <= x, z >= 0
            x, y = np.abs(x), np.abs(y)
            x, y, z = np.maximum(x, y), np.minimum(x, y), np.abs(z)
        elif self.lattice == 'orthorhombic':
            x, y, z = np.abs(x), np.abs(y), np.abs(z)
        elif self.lattice == 'monoclinic':  # 2-fold axis along y: (-x, y, -z), mirror (x, -y, z)
            z = np.where(x < 0, -z, z)
            x, y = np.abs(x), np.abs(y)
        elif self.lattice == 'triclinic':
            x, y, z = np.where(z < 0, [-x, -y, -z], [x, y, z])
        elif self.lattice in ('hexagonal', 'trigonal'):
            # inversion into the upper hemisphere, then rotations and mirrors fold the azimuth into the wedge
            # [offset, offset + period/2]: hexagonal 0-30°, trigonal 30-90°
            period, offset = (np.pi/3, 0.) if self.lattice == 'hexagonal' else (2*np.pi/3, np.pi/6)
            x, y, z = np.where(z < 0, [-x, -y, -z], [x, y, z])
            azimuth = np.mod(np.arctan2(y, x) - offset, period)
            azimuth = np.minimum(azimuth, period - azimuth) + offset
            radius = np.hypot(x, y)
            x, y = radius*np.cos(azimuth), radius*np.sin(azimuth)
        return np.array([x, y, z])


    def xyToHKL(self, inPlane: Any) -> np.ndarray | None:
        """Convert cubic stereographic-plane coordinates to HKL vectors.

        Args:
           inPlane: coordinates (2,) or (2, n) in the stereographic plane

        Returns:
           HKL vectors (3,) or (3, n); None if the lattice is not cubic
        """
        if self.lattice != 'cubic':
            print('ERROR: only implemented for cubic lattice')
            return None
        inPlane = np.asarray(inPlane, dtype=float)
        l = 2.0 / (inPlane[0]*inPlane[0] + inPlane[1]*inPlane[1] + 1.)
        return np.concatenate((l*inPlane, np.asarray(l-1.0)[None]))


    def __repr__(self) -> str:
        """Return the lattice name.

        Returns:
           lattice name
        """
        return str(self.lattice)


    def __eq__(self, other: object) -> bool:
        """Compare lattices.

        Args:
           other: object to compare with

        Returns:
           True if both lattices are equal
        """
        if not isinstance(other, Symmetry):
            return NotImplemented
        return self.lattice == other.lattice


def _hexagonalCell(a: float, c: float) -> np.ndarray:
    """Return the eighteen centered edges of a regular hexagonal prism.

    Args:
       a: lattice constant a
       c: lattice constant c

    Returns:
       edges (18, 6)
    """
    angles = np.arange(6) * np.pi / 3.0
    basal = np.column_stack((a * np.cos(angles), a * np.sin(angles)))
    lower = np.column_stack((basal, np.full(6, -c / 2.0)))
    upper = np.column_stack((basal, np.full(6, c / 2.0)))
    vertices = np.vstack((lower, upper))
    connections = ([(index, (index + 1) % 6) for index in range(6)] +
                   [(index + 6, (index + 1) % 6 + 6) for index in range(6)] +
                   [(index, index + 6) for index in range(6)])
    return np.array([np.concatenate((vertices[start], vertices[end])) for start, end in connections])


def _parallelepipedCell(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Return centered edges of a cell defined by lengths and angles in degrees.

    Args:
       a: lattice constant a
       b: lattice constant b
       c: lattice constant c
       alpha: angle between b and c in degrees
       beta: angle between a and c in degrees
       gamma: angle between a and b in degrees

    Returns:
       edges (12, 6)
    """
    alpha, beta, gamma = np.deg2rad((alpha, beta, gamma))
    cosAlpha, cosBeta, cosGamma = np.cos((alpha, beta, gamma))
    sinGamma = np.sin(gamma)
    cY = c * (cosAlpha - cosBeta * cosGamma) / sinGamma
    cZSquared = c*c - (c*cosBeta)**2 - cY*cY
    if cZSquared <= 0.0:
        raise ValueError('lattice angles do not define a valid unit cell')
    vectors = np.array([
        [a, 0.0, 0.0],
        [b * cosGamma, b * sinGamma, 0.0],
        [c * cosBeta, cY, np.sqrt(cZSquared)]])
    vertices = np.array([(sx * vectors[0] + sy * vectors[1] + sz * vectors[2]) / 2.0
                         for sx in (-1.0, 1.0) for sy in (-1.0, 1.0) for sz in (-1.0, 1.0)])
    connections = ((0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
                   (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7))
    return np.array([np.concatenate((vertices[start], vertices[end])) for start, end in connections])


def showOrSave(fileName: str = '', show: bool = True) -> None:
    """Save the current figure to a file, or else show it.

    Args:
       fileName: if given, save to this file and close the figure
       show: show the figure if no fileName is given
    """
    if fileName:
        plt.savefig(fileName, dpi=150, bbox_inches='tight')
        plt.close()
    elif show:
        plt.show()
