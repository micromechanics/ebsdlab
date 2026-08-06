#!/usr/bin/python3
"""Tests for material symmetry operations and visualization geometry."""

import numpy as np
import pytest
from ebsdlab.symmetry import GROUPS, Symmetry


def test_group_configuration_is_complete():
    """Every lattice is fully described by the single group registry."""
    assert set(GROUPS) == {
        'cubic', 'hexagonal', 'tetragonal', 'orthorhombic'
    }
    for config in GROUPS.values():
        assert set(config) == {'rotation_group', 'sst_bases', 'cell'}
        assert set(config['sst_bases']) == {'improper', 'proper'}
        assert set(config['cell']) == {
            'geometry', 'default_ratio', 'equal_to_a'
        }
        assert config['cell']['geometry'] in {'orthogonal', 'hexagonal'}
        for basis in config['sst_bases'].values():
            assert basis.shape == (3, 3)


def test_proper_and_improper_standard_stereographic_triangles():
    """Improper symmetry folds opposite poles; proper symmetry needs two SSTs."""
    symmetry = Symmetry('cubic')

    assert symmetry.inSST([0.0, 0.0, -1.0])
    assert not symmetry.inSST([0.0, 0.0, -1.0], proper=True)
    assert not symmetry.inSST([0.0, 1.0, 1.0])
    assert symmetry.inSST([0.0, 1.0, 1.0], proper=True)


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


@pytest.mark.parametrize(
    ('lattice', 'shape', 'lengths'),
    [
        ('cubic', (12, 6), [1.0] * 12),
        ('tetragonal', (12, 6), [1.0] * 8 + [1.5] * 4),
        ('hexagonal', (18, 6), [1.0] * 12 + [1.5] * 6),
        ('orthorhombic', (12, 6), [1.0] * 4 + [1.25] * 4 + [1.5] * 4),
    ],
)
def test_default_unit_cells(lattice, shape, lengths):
    """Every supported lattice has centered, finite, illustrative geometry."""
    cell = Symmetry(lattice).unitCell()
    edge_lengths = np.linalg.norm(cell[:, 3:] - cell[:, :3], axis=1)

    assert cell.shape == shape
    assert np.all(np.isfinite(cell))
    np.testing.assert_allclose(
        np.mean(cell.reshape(-1, 3), axis=0), 0.0, atol=1e-15
    )
    np.testing.assert_allclose(np.sort(edge_lengths), np.sort(lengths), atol=1e-15)


@pytest.mark.parametrize(
    ('lattice', 'parameters', 'lengths'),
    [
        ('cubic', {'a': 2.0}, [2.0] * 12),
        ('tetragonal', {'a': 2.0, 'c': 3.0}, [2.0] * 8 + [3.0] * 4),
        ('hexagonal', {'a': 2.0, 'c': 3.0}, [2.0] * 12 + [3.0] * 6),
        ('orthorhombic', {'a': 2.0, 'b': 3.0, 'c': 4.0},
         [2.0] * 4 + [3.0] * 4 + [4.0] * 4),
    ],
)
def test_unit_cells_accept_lattice_constants(lattice, parameters, lengths):
    """Explicit lattice constants determine the generated edge lengths."""
    cell = Symmetry(lattice).unitCell(**parameters)
    edge_lengths = np.linalg.norm(cell[:, 3:] - cell[:, :3], axis=1)

    np.testing.assert_allclose(np.sort(edge_lengths), np.sort(lengths), atol=1e-14)


@pytest.mark.parametrize(
    ('lattice', 'parameters'),
    [
        ('cubic', {'a': 0.0}),
        ('cubic', {'c': 2.0}),
        ('tetragonal', {'b': 2.0}),
        ('hexagonal', {'b': 2.0}),
        ('orthorhombic', {'c': np.inf}),
    ],
)
def test_unit_cells_reject_invalid_lattice_constants(lattice, parameters):
    """Dimensions must be physical and respect constrained equal axes."""
    with pytest.raises(ValueError):
        Symmetry(lattice).unitCell(**parameters)
