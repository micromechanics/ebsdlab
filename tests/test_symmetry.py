#!/usr/bin/python3
"""TEST Symmetry class """
import numpy as np
import pytest
from ebsdlab.symmetry import Symmetry


@pytest.mark.mpl_image_compare
def test_symmetry():
    s = Symmetry('cubic')
    point1 = np.array((2., 1., 3.))
    point1 /= np.linalg.norm(point1)
    inside, color = s.inSST(point1, color=True)
    assert inside, 'point is inside the SST'
    assert np.linalg.norm(color-np.array([0.75983569, 0.903602, 1.])
                          ) < 0.0001, 'Color should be [0.75983569,0.903602,1.]'

    point2 = np.array((1., 2., 3.))
    point2 /= np.linalg.norm(point2)
    inside, color = s.inSST(point2, color=True)
    assert not inside, 'point is outside the SST'
    assert np.linalg.norm(color) < 0.0001, 'Color should be black'

    fig = s.standardTriangle()
    return fig
