"""Crystal symmetry operations and visualization helpers."""

import math
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.spatial.transform import Rotation

from ._rotation import as_rodrigues

# standard stereographic triangle (SST): Smallest region representing symmetry-equivalent crystal directions in an inverse pole figure (IPF).
# The improper basis uses inversion/reflection equivalence and folds opposite directions into one triangle.
# The proper basis describes the adjoining triangle needed when equivalence is restricted to handedness-preserving rotations.
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
        },
    },
}


def _orthogonal_cell(a: float, b: float, c: float) -> np.ndarray:
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


def _hexagonal_cell(a: float, c: float) -> np.ndarray:
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



class Symmetry:
    """Material symmetry identified by a human-readable lattice name."""

    def __init__(self, symmetry: str | None = None) -> None:
        self.lattice: str | None
        if isinstance(symmetry, str) and symmetry.lower() in GROUPS:
            self.lattice = symmetry.lower()
        else:
            self.lattice = None

    def __copy__(self) -> 'Symmetry':
        return self.__class__(self.lattice)
    copy = __copy__

    def __repr__(self) -> str:
        return str(self.lattice)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Symmetry):
            return NotImplemented
        return self.lattice == other.lattice

    def symmetryQuats(self, who: Any = None) -> Rotation:
        """Return the proper symmetry rotations, with identity at index zero."""
        # Generate all symmetry-equivalent rotation operations
        config = GROUPS.get(self.lattice) if self.lattice is not None else None
        operations = (Rotation.create_group(config['rotation_group']) if config is not None else Rotation.identity(1))
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

    def unitCell(self, *, a: float = 1.0, b: float | None = None, c: float | None = None) -> np.ndarray | list[list[None]]:
        """Return centered unit-cell edges as ``[x1, y1, z1, x2, y2, z2]``.

        The defaults are illustrative proportions for plotting. Supply positive
        lattice constants for a material-specific cell.
        """
        if self.lattice not in GROUPS:
            print('Unit cell not implemented')
            return [[None]]

        cellConfig = GROUPS[self.lattice]['cell']
        _, default_b, default_c = cellConfig['default_ratio']
        b = a * default_b if b is None else b
        c = a * default_c if c is None else c
        dimensions = np.asarray((a, b, c), dtype=float)
        if not np.all(np.isfinite(dimensions)) or np.any(dimensions <= 0.0):
            raise ValueError('lattice constants must be finite and positive')
        axisValues = {'b': b, 'c': c}
        for axis in cellConfig['equal_to_a']:
            if not np.isclose(a, axisValues[axis]):
                raise ValueError(f'{self.lattice} requires a = {axis}')
        if cellConfig['geometry'] == 'hexagonal':
            return _hexagonal_cell(a, c)
        return _orthogonal_cell(a, b, c)


    def equivalentQuaternions(self, quaternion: Rotation, who: Any = None) -> list[Rotation]:
        """Return rotations equivalent to ``quaternion`` under this symmetry."""
        return [quaternion*q for q in self.symmetryQuats(who)]


    def inFZ(self, R: Rotation | np.ndarray) -> bool:
        """Return whether a Rodrigues vector lies in the fundamental zone."""
        if isinstance(R, Rotation):
            R = as_rodrigues(R)
        # fundamental zone in Rodrigues space is point symmetric around origin
        R = abs(R)
        if self.lattice == 'cubic':
            limit = math.sqrt(2.0) - 1.0
            return bool(
                limit >= R[0] and limit >= R[1] and limit >= R[2]
                and 1.0 >= R[0] + R[1] + R[2]
            )
        if self.lattice == 'hexagonal':
            return bool(
                1.0 >= R[0] and 1.0 >= R[1] and 1.0 >= R[2]
                and 2.0 >= math.sqrt(3.0)*R[0] + R[1]
                and 2.0 >= math.sqrt(3.0)*R[1] + R[0]
                and 2.0 >= math.sqrt(3.0) + R[2]
            )
        if self.lattice == 'tetragonal':
            return bool(
                1.0 >= R[0] and 1.0 >= R[1]
                and math.sqrt(2.0) >= R[0] + R[1]
                and math.sqrt(2.0) >= R[2] + 1.0
            )
        if self.lattice == 'orthorhombic':
            return bool(1.0 >= R[0] and 1.0 >= R[1] and 1.0 >= R[2])
        return True


    def inDisorientationSST(self, R: Rotation | np.ndarray) -> bool:
        """Return whether a misorientation lies in the standard triangle.

        The criteria follow Heinz and Neumann, Acta Cryst. A47 (1991), 780-789.
        """
        if isinstance(R, Rotation):
            R = as_rodrigues(R)
        if self.lattice == 'cubic':
            return bool(R[0] >= R[1] and R[1] >= R[2] and R[2] >= 0.0)
        if self.lattice == 'hexagonal':
            return bool(
                R[0] >= math.sqrt(3.0)*R[1]
                and R[1] >= 0.0 and R[2] >= 0.0
            )
        if self.lattice == 'tetragonal':
            return bool(R[0] >= R[1] and R[1] >= 0.0 and R[2] >= 0.0)
        if self.lattice == 'orthorhombic':
            return bool(R[0] >= 0.0 and R[1] >= 0.0 and R[2] >= 0.0)
        return True


    def inSST(self, vector: Any, proper: bool = False, color: bool = False) -> Any:
        """Test vectors against the standard triangle and optionally color them.

        ``proper`` also considers the neighboring proper triangle. With
        ``color=True``, return the membership flags and IPF colors as RGB values.
        """
        config = GROUPS.get(self.lattice) if self.lattice is not None else None
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
                inSST = np.all(theComponents >= 0.0, axis=0)
                # ... and proper SST
                if not np.any(inSST):
                    theComponents = np.dot(basis['proper'], v)
            else:
                # z component projects identical for positive and negative values
                v[2] = abs(v[2])
                theComponents = np.dot(basis['improper'], v)
            inSST = np.all(theComponents >= 0.0, axis=0)
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


    def standardTriangle(self, fileName: str | None = None, show: bool = True, stepSize: float = 0.1) -> Any:
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
