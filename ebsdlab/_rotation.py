"""Convention adapters around :mod:`scipy.spatial.transform`."""

import warnings
import numpy as np
from scipy.spatial.transform import Rotation


def asRodrigues(rotation: Rotation) -> np.ndarray:
    """Return Rodrigues vectors ``axis * tan(angle / 2)``.

    Args:
       rotation: one or several rotations

    Returns:
       Rodrigues vectors of shape (..., 3); infinite for 180° rotations
    """
    quaternion = rotation.as_quat(scalar_first=True)
    quaternion = np.where(quaternion[..., :1] < 0.0, -quaternion, quaternion)
    with np.errstate(divide='ignore', invalid='ignore'):
        return quaternion[..., 1:] / quaternion[..., :1]


def asBungeEulers(rotation: Rotation, degrees: bool = False, standardRange: bool = False) -> np.ndarray:
    """Return intrinsic active ZXZ Euler angles in the Bunge convention.

    Args:
       rotation: one or several rotations
       degrees: return degrees instead of radians
       standardRange: map phi1 and phi2 into [0, 2pi)

    Returns:
       Euler angles (phi1, Phi, phi2) of shape (..., 3)
    """
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Gimbal lock detected')
        angles = rotation.as_euler('ZXZ', degrees=False)
    if standardRange:
        angles = np.array(angles, copy=True)
        angles[..., [0, 2]] %= 2.0 * np.pi
        angles[np.isclose(angles, 2.0 * np.pi, atol=1e-12)] = 0.0
    return np.degrees(angles) if degrees else angles


# For large arrays, the quaternion arithmetic below is several times faster than the same scipy Rotation method.
def fromBungeEulers(eulers: np.ndarray) -> Rotation:
    """Return rotations from Bunge Euler angles, as ``Rotation.from_euler('ZXZ', eulers)``.

    Args:
       eulers: Euler angles (phi1, Phi, phi2) in radians, shape (3,) or (n, 3)

    Returns:
       rotations
    """
    phi1, phi, phi2 = np.moveaxis(np.asarray(eulers, dtype=float), -1, 0)
    halfSum, halfDiff = (phi1+phi2)/2, (phi1-phi2)/2
    return Rotation.from_quat(np.stack([np.sin(phi/2)*np.cos(halfDiff), np.sin(phi/2)*np.sin(halfDiff),
                                        np.cos(phi/2)*np.sin(halfSum), np.cos(phi/2)*np.cos(halfSum)], axis=-1))


def multiply(first: Rotation, second: Rotation) -> Rotation:
    """Return the product ``first * second``: first apply second, then first.

    Args:
       first: one or several rotations
       second: one or several rotations, broadcast with first

    Returns:
       rotations
    """
    x1, y1, z1, w1 = np.moveaxis(first.as_quat(), -1, 0)
    x2, y2, z2, w2 = np.moveaxis(second.as_quat(), -1, 0)
    return Rotation.from_quat(np.stack([w1*x2 + x1*w2 + y1*z2 - z1*y2, w1*y2 - x1*z2 + y1*w2 + z1*x2,
                                        w1*z2 + x1*y2 - y1*x2 + z1*w2, w1*w2 - x1*x2 - y1*y2 - z1*z2], axis=-1))
