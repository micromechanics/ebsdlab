"""Small convention adapters around :mod:`scipy.spatial.transform`."""

import warnings
import numpy as np
from scipy.spatial.transform import Rotation


def asRodrigues(rotation: Rotation) -> np.ndarray:
    """Return Rodrigues vectors ``axis * tan(angle / 2)``."""
    quaternion = rotation.as_quat(scalar_first=True)
    quaternion = np.where(quaternion[..., :1] < 0.0, -quaternion, quaternion)
    with np.errstate(divide="ignore", invalid="ignore"):
        return quaternion[..., 1:] / quaternion[..., :1]


def asBungeEulers(
    rotation: Rotation, degrees: bool = False, standardRange: bool = False
) -> np.ndarray:
    """Return intrinsic active ZXZ Euler angles in the Bunge convention."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Gimbal lock detected")
        angles = rotation.as_euler("ZXZ", degrees=False)
    if standardRange:
        angles = np.array(angles, copy=True)
        angles[..., [0, 2]] %= 2.0 * np.pi
        angles[np.isclose(angles, 2.0 * np.pi, atol=1e-12)] = 0.0
    return np.degrees(angles) if degrees else angles
