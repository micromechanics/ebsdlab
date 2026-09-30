"""Crystal symmetry operations and visualization helpers."""

import math
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.spatial.transform import Rotation

from ._rotation import asRodrigues

# Directions on an SST edge give components of ~1e-17 instead of exactly 0; count them as inside
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
            # A 60 degree azimuthal wedge in the upper/lower hemisphere.
            'improper': np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0],
                                  [np.sqrt(3.0), -1.0, 0.0]]),
            'proper': np.array([[0.0, 0.0, -1.0], [0.0, 1.0, 0.0],
                                [np.sqrt(3.0), -1.0, 0.0]]),
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

# Rhombohedral is the lattice name customarily used for the trigonal crystal
# system.  Keep one canonical registry entry so both names behave identically.
LATTICE_ALIASES = {'rhombohedral': 'trigonal'}


def _orthogonalCell(a: float, b: float, c: float) -> np.ndarray:
    """Return the twelve centered edges of an orthogonal unit cell."""
    vertices = np.array([
        [-a, -b, -c], [-a, -b, c], [-a, b, -c], [-a, b, c],
        [a, -b, -c], [a, -b, c], [a, b, -c], [a, b, c],
    ]) / 2.0
    connections = (
        (0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
        (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7),
    )
    return np.array([np.concatenate((vertices[start], vertices[end])) for start, end in connections])


def _hexagonalCell(a: float, c: float) -> np.ndarray:
    """Return the eighteen centered edges of a regular hexagonal prism."""
    angles = np.arange(6) * np.pi / 3.0
    basal = np.column_stack((a * np.cos(angles), a * np.sin(angles)))
    lower = np.column_stack((basal, np.full(6, -c / 2.0)))
    upper = np.column_stack((basal, np.full(6, c / 2.0)))
    vertices = np.vstack((lower, upper))
    connections = (
        [(index, (index + 1) % 6) for index in range(6)]
        + [(index + 6, (index + 1) % 6 + 6) for index in range(6)]
        + [(index, index + 6) for index in range(6)]
    )
    return np.array([np.concatenate((vertices[start], vertices[end])) for start, end in connections])


def _parallelepipedCell(a: float, b: float, c: float,
                         alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Return centered edges of a cell defined by lengths and angles in degrees."""
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
        [c * cosBeta, cY, np.sqrt(cZSquared)],
    ])
    vertices = np.array([
        (sx * vectors[0] + sy * vectors[1] + sz * vectors[2]) / 2.0
        for sx in (-1.0, 1.0) for sy in (-1.0, 1.0) for sz in (-1.0, 1.0)
    ])
    connections = ((0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
                   (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7))
    return np.array([np.concatenate((vertices[start], vertices[end])) for start, end in connections])



class Symmetry:
    """Material symmetry identified by a human-readable lattice name."""

    def __init__(self, symmetry: str = '') -> None:
        """Create a symmetry from a lattice name such as 'cubic'; '' means not identified."""
        self.lattice = ''
        if symmetry:
            lattice = LATTICE_ALIASES.get(symmetry.lower(), symmetry.lower())
            if lattice not in GROUPS:
                raise ValueError(f'Unknown symmetry {symmetry!r}. Supported: {", ".join(GROUPS)}')
            self.lattice = lattice

    def __copy__(self) -> 'Symmetry':
        return self.__class__(self.lattice)

    def __repr__(self) -> str:
        return str(self.lattice)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Symmetry):
            return NotImplemented
        return self.lattice == other.lattice

    def symmetryQuats(self, who: Any = None) -> Rotation:
        """Return the proper symmetry rotations, with identity at index zero."""
        # Generate all symmetry-equivalent rotation operations
        config = GROUPS.get(self.lattice)
        operations = (Rotation.create_group(
            config['rotation_group'], axis=config.get('rotation_axis', 'Z')
        ) if config is not None else Rotation.identity(1))
        # SciPy's cubic group does not start with identity, but callers rely on
        # operations[0] being the rotation that leaves vectors unchanged.
        identity = int(np.argmin(operations.magnitude()))
        if identity != 0:
            order = np.concatenate(([identity], np.delete(np.arange(len(operations)), identity)))
            operations = operations[order]
        # Select a subset of those equivalent rotations
        if who is not None and np.size(who) > 0:
            operations = operations[np.atleast_1d(who)]
        return operations

    def unitCell(self, *, a: float = 1.0, b: float | None = None,
                 c: float | None = None, alpha: float | None = None,
                 beta: float | None = None, gamma: float | None = None) -> np.ndarray | list[list[None]]:
        """Return centered unit-cell edges as ``[x1, y1, z1, x2, y2, z2]``.

        The defaults are illustrative proportions for plotting. Supply positive
        lattice constants for a material-specific cell.
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
        angles = np.asarray(tuple(
            default if value is None else value
            for value, default in zip((alpha, beta, gamma), defaultAngles)
        ), dtype=float)
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
        if cellConfig['geometry'] in {'monoclinic', 'triclinic', 'rhombohedral'}:
            return _parallelepipedCell(a, b, c, *angles)
        return _orthogonalCell(a, b, c)


    def equivalentQuaternions(self, quaternion: Rotation, who: Any = None) -> list[Rotation]:
        """Return rotations equivalent to ``quaternion`` under this symmetry."""
        return [quaternion*q for q in self.symmetryQuats(who)]


    def inFZ(self, rotation: Rotation) -> bool:
        """Return whether a rotation lies in the fundamental zone (FZ)"""
        rawRodrigues = np.asarray(asRodrigues(rotation), dtype=float)
        # fundamental zone in Rodrigues space is point symmetric around origin
        rodrigues = abs(rawRodrigues)
        if self.lattice == 'cubic':
            limit = math.sqrt(2.0) - 1.0
            return bool(
                limit >= rodrigues[0] and limit >= rodrigues[1] and limit >= rodrigues[2]
                and 1.0 >= rodrigues[0] + rodrigues[1] + rodrigues[2]
            )
        if self.lattice == 'hexagonal':
            return bool(
                1.0 >= rodrigues[0] and 1.0 >= rodrigues[1] and 1.0 >= rodrigues[2]
                and 2.0 >= math.sqrt(3.0)*rodrigues[0] + rodrigues[1]
                and 2.0 >= math.sqrt(3.0)*rodrigues[1] + rodrigues[0]
                and 2.0 >= math.sqrt(3.0) + rodrigues[2]
            )
        if self.lattice == 'tetragonal':
            return bool(
                1.0 >= rodrigues[0] and 1.0 >= rodrigues[1]
                and math.sqrt(2.0) >= rodrigues[0] + rodrigues[1]
                and math.sqrt(2.0) >= rodrigues[2] + 1.0
            )
        if self.lattice == 'orthorhombic':
            return bool(1.0 >= rodrigues[0] and 1.0 >= rodrigues[1] and 1.0 >= rodrigues[2])
        if self.lattice in {'monoclinic', 'triclinic', 'trigonal'}:
            # The Voronoi region of identity is the fundamental zone for
            # these lower-symmetry proper rotation groups.
            magnitude = np.linalg.norm(rawRodrigues)
            rotation = (Rotation.from_rotvec(
                2.0 * math.atan(magnitude) * rawRodrigues / magnitude
            ) if magnitude else Rotation.identity())
            magnitudes = (rotation * self.symmetryQuats()).magnitude()
            return bool(magnitudes[0] <= np.min(magnitudes) + 1e-12)
        return True


    def inDisorientationSST(self, rotation: Rotation) -> bool:
        """Return whether a misorientation lies in the standard triangle.

        The criteria follow Heinz and Neumann, Acta Cryst. A47 (1991), 780-789.
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
                inSST = inSST | oppositeSST
            inSST = np.all(theComponents >= -SST_TOLERANCE, axis=0)
        # have to return color array
        if color:
            if np.any(inSST):
                theComponentsNorm = (
                    theComponents / np.linalg.norm(theComponents, axis=0)
                )
                # smoothen color ramps
                rgb = np.power(np.abs(theComponentsNorm), 0.5)
                if rgb.ndim > 1:
                    rgb[:, ~inSST] = 0.0
                # limit to maximum intensity
                rgb[rgb > 1.0] = 1.0
                if rgb.ndim > 1:
                    # normalize to (HS)V = 1
                    rgb[:, inSST] /= np.max(rgb[:, inSST], axis=0)
                else:
                    rgb /= np.max(rgb)
            else:
                rgb = np.zeros(3, 'd')
            return (inSST, rgb)
        return inSST


    def xyToHKL(self, inPlane: Any) -> np.ndarray | None:
        """Convert cubic stereographic-plane coordinates to HKL vectors."""
        if self.lattice != 'cubic':
            print('ERROR: only implemented for cubic lattice')
            return None
        inPlane = np.asarray(inPlane, dtype=float)
        l = 2.0 / (inPlane[0]*inPlane[0] + inPlane[1]*inPlane[1] + 1.)
        if inPlane.ndim == 1:
            hkl = np.zeros(3, dtype=float)
            hkl[:2] = l*inPlane
            hkl[2] = l-1.0
        else:
            hkl = np.zeros((3, inPlane.shape[1]), dtype=float)
            hkl[:2, :] = l*inPlane
            hkl[2, :] = l-1.0
        return hkl


    def standardTriangle(self, fileName: str = '', show: bool = True, stepSize: float = 0.1) -> Any:
        """Plot the colored cubic standard stereographic triangle."""
        if self.lattice != 'cubic':
            print('ERROR: only implemented for cubic lattice')
            return None
        # create border
        borderPoints = [[0.0, 0.0]]
        xTemp = []
        yTemp = []
        for a in range(16):
            borderPoints.append([
                math.cos(a/180.0*math.pi)*math.sqrt(2.0) - 1,
                math.sin(a/180.0*math.pi)*math.sqrt(2.0),
            ])
            xTemp.append(math.cos(a/180.0*math.pi)*math.sqrt(2.0)-1)
            yTemp.append(math.sin(a/180.0*math.pi)*math.sqrt(2))
        borderPoints.append([0.0, 0.0])
        border = np.array(borderPoints)
        func = interp1d(yTemp, xTemp)  # function yTemp=func(xTemp)
        # create colored background
        xyPoints = []
        for iY in np.arange(0, 0.366025403784, stepSize):
            for iX in np.arange(iY, 0.41421, stepSize):
                if func(iY) <= iX:
                    continue
                xyPoints.append([iX, iY])
        xy = np.array(xyPoints)

        hkl = self.xyToHKL(xy.T)
        assert hkl is not None
        _, rgb = self.inSST(hkl, color=True, proper=False)
        colors = []
        for i in range(rgb.shape[1]):
            values = (rgb[:, i]*255).astype(int)
            string = f'#{values[0]:02x}{values[1]:02x}{values[2]:02x}'
            colors.append(string)
        # plotting: background, border, labels
        fig, ax = plt.subplots()
        ax.scatter(xy[:, 0], xy[:, 1], c=colors, s=15000. *
                   stepSize, linewidths=0)  # , alpha=0.05)
        ax.plot(border[:, 0], border[:, 1], '-k', linewidth=3)  # , alpha=0.5)
        plt.rcParams['font.size'] = 18.
        ax.text(0, 0, '[100]', horizontalalignment='right')  # , zorder=40
        ax.text(0.42, 0, '[110]')
        ax.text(0.37, 0.37, '[111]')
        ax.axis('off')
        ax.axis('equal')
        if fileName:
            plt.savefig(fileName)
        if show:
            plt.show()
        return fig
