"""Lock the rotation and crystallographic conventions used by ebsdlab."""

import numpy as np
from scipy.spatial.transform import Rotation

from ebsdlab.ebsd import EBSD
from ebsdlab.orientation import Orientation
from ebsdlab.symmetry import Symmetry


ATOL = 1e-12


def _rotation_x(angle):
    """Return an active right-handed rotation about x."""
    cosine = np.cos(angle)
    sine = np.sin(angle)
    return np.array(
        [[1.0, 0.0, 0.0], [0.0, cosine, -sine], [0.0, sine, cosine]]
    )


def _rotation_z(angle):
    """Return an active right-handed rotation about z."""
    cosine = np.cos(angle)
    sine = np.sin(angle)
    return np.array(
        [[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]]
    )


def test_quaternion_is_scalar_first_and_rotates_actively():
    """A positive z rotation is stored as wxyz and maps +x onto +y."""
    quaternion = np.array([np.sqrt(0.5), 0.0, 0.0, np.sqrt(0.5)])
    rotation = Rotation.from_quat(quaternion, scalar_first=True)

    np.testing.assert_allclose(
        rotation.as_quat(scalar_first=True),
        quaternion,
        atol=ATOL,
    )
    np.testing.assert_allclose(
        rotation.apply([1.0, 0.0, 0.0]),
        [0.0, 1.0, 0.0],
        atol=ATOL,
    )


def test_bunge_eulers_are_intrinsic_active_zxz():
    """Bunge angles mean the active intrinsic ZXZ convention."""
    phi1, phi, phi2 = np.deg2rad([10.0, 20.0, 30.0])
    orientation = Orientation(eulers=np.array([phi1, phi, phi2]))
    expected = _rotation_z(phi1) @ _rotation_x(phi) @ _rotation_z(phi2)

    np.testing.assert_allclose(orientation.asMatrix(), expected, atol=ATOL)
    np.testing.assert_allclose(
        orientation.asEulers(standardRange=True),
        [phi1, phi, phi2],
        atol=ATOL,
    )


def test_quaternion_composition_applies_right_operand_first():
    """For qz*qx, qx acts on a vector before qz."""
    rotate_x = Rotation.from_rotvec(np.pi / 2.0 * np.array([1.0, 0.0, 0.0]))
    rotate_z = Rotation.from_rotvec(np.pi / 2.0 * np.array([0.0, 0.0, 1.0]))
    vector = np.array([0.0, 1.0, 0.0])

    composed = (rotate_z * rotate_x).apply(vector)
    sequential = rotate_z.apply(rotate_x.apply(vector))

    np.testing.assert_allclose(composed, sequential, atol=ATOL)
    np.testing.assert_allclose(composed, [0.0, 0.0, 1.0], atol=ATOL)


def test_inverse_pole_uses_the_inverse_orientation():
    """A sample direction is mapped back into the crystal frame."""
    orientation = Orientation(
        eulers=np.deg2rad(np.array([0.0, 90.0, 0.0])), symmetry=None
    )

    pole, symmetry_index = orientation.inversePole(
        np.array([0.0, 0.0, 1.0]), sst=False
    )

    np.testing.assert_allclose(pole, [0.0, 1.0, 0.0], atol=ATOL)
    assert symmetry_index == 0


def test_cubic_disorientation_is_symmetry_reduced():
    """A 100 degree z rotation is 10 degrees from a cubic equivalent."""
    identity = Orientation(eulers=np.zeros(3), symmetry="cubic")
    rotated = Orientation(
        eulers=np.deg2rad(np.array([100.0, 0.0, 0.0])), symmetry="cubic"
    )

    disorientation, _, _, _ = identity.disorientation(rotated)
    angle = np.degrees(disorientation.quaternion.magnitude())

    np.testing.assert_allclose(angle, 10.0, atol=1e-10)


def test_symmetry_groups_start_with_identity():
    """Algorithms and public indexing use operation zero as the identity."""
    expected_sizes = {
        "cubic": 24,
        "hexagonal": 12,
        "tetragonal": 8,
        "orthorhombic": 4,
        "monoclinic": 2,
        "triclinic": 1,
        "trigonal": 6,
    }

    for lattice, size in expected_sizes.items():
        operations = Symmetry(lattice).symmetryQuats()
        assert len(operations) == size
        np.testing.assert_allclose(operations[0].as_matrix(), np.eye(3), atol=ATOL)


def test_documented_cubic_orientation_average():
    """Lock the symmetry-aware average used by the documentation examples."""
    orientations = [
        Orientation(eulers=np.deg2rad(angles), symmetry="cubic")
        for angles in ([0.0, 45.0, 0.0], [0.0, 0.0, 0.0], [0.0, 15.0, 0.0])
    ]

    average = Orientation.average(orientations)

    np.testing.assert_allclose(
        average.asEulers(degrees=True, standardRange=True),
        [0.0, 19.86780516, 0.0],
        atol=1e-8,
    )


def test_kam_is_mean_symmetry_reduced_misorientation():
    """A point surrounded by six ten-degree neighbors has a KAM of ten."""
    ebsd = EBSD.__new__(EBSD)
    ebsd.x = np.arange(7, dtype=float)
    ebsd.ci = np.ones(7)
    ebsd.sym = [Symmetry("cubic")]

    eulers = np.zeros((3, 7))
    eulers[1, 1:] = np.deg2rad(10.0)
    ebsd.quaternions = Rotation.from_euler("ZXZ", eulers.T)

    neighbors = np.repeat(np.arange(7)[:, np.newaxis], 6, axis=1)
    neighbors[0] = np.arange(1, 7)
    ebsd.neighbors = lambda: neighbors

    ebsd.calcKAM()

    np.testing.assert_allclose(ebsd.kam[0], 10.0, atol=1e-10)
    np.testing.assert_allclose(ebsd.kam[1:], 0.0, atol=ATOL)


def test_cubic_ipf_colors_at_high_symmetry_directions():
    """The cubic [001], [101], and [111] corners are red, green, blue."""
    red = Orientation(eulers=np.zeros(3), symmetry="cubic")
    green = Orientation(
        eulers=np.array([0.0, np.pi / 4.0, 0.0]), symmetry="cubic"
    )
    blue = Orientation(
        eulers=np.array([0.0, np.arccos(1.0 / np.sqrt(3.0)), np.pi / 4.0]),
        symmetry="cubic",
    )
    sample_normal = [0.0, 0.0, 1.0]

    np.testing.assert_allclose(red.ipfColor(sample_normal), [1.0, 0.0, 0.0], atol=ATOL)
    np.testing.assert_allclose(
        green.ipfColor(sample_normal), [0.0, 1.0, 0.0], atol=2e-8
    )
    np.testing.assert_allclose(
        blue.ipfColor(sample_normal), [0.0, 0.0, 1.0], atol=2e-8
    )
