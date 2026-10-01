#!/usr/bin/python3
"""TEST EBSD class """
from pathlib import Path
import numpy as np
import pytest
from ebsdlab.ebsd import EBSD, SUPPORTED_SUFFIXES


DATA_DIR = Path(__file__).parent/'DataFiles'
ANG_DATA = '''# Symmetry 62
# OPERATOR
0 0 0 0 0 1 1 1 1 1
0 0 0 1 0 1 1 1 1 1
0 0 0 0 1 1 1 1 1 1
0 0 0 1 1 1 1 1 1 1
'''


@pytest.mark.parametrize('data_file',
    sorted(path for path in DATA_DIR.rglob('*') if path.is_file() and path.suffix.lower() in SUPPORTED_SUFFIXES),
    ids=lambda path: path.name)
def test_all_supported_data_files_load(data_file):
    """Every supported fixture shipped with the tests can be opened."""
    ebsd = EBSD(str(data_file))
    assert len(ebsd.x) > 0
    assert len(ebsd.quaternions) == len(ebsd.x)


def test_osc_ipf_with_manually_supplied_symmetry():
    """OSC maps can use an explicitly supplied phase symmetry for IPF colors."""
    ebsd = EBSD(str(DATA_DIR/'EBSD.osc'), symmetry='cubic')
    ebsd.plotIPF(show=False)

    assert repr(ebsd.sym[0]) == 'cubic'
    assert np.any(np.asarray(ebsd.image))


def test_ang_with_manual_non_cubic_symmetry_loads(tmp_path):
    """A manually supplied symmetry permits ANG data with unknown metadata."""
    data_file = tmp_path/'hexagonal.ang'
    data_file.write_text(ANG_DATA)
    ebsd = EBSD(data_file, symmetry='hexagonal')
    assert len(ebsd.x) == 4
    assert repr(ebsd.sym[0]) == 'hexagonal'


def test_file_extensions_are_case_insensitive_and_invalid_ones_raise(tmp_path):
    """Bad user input must not terminate the embedding Python process."""
    data_file = tmp_path/'map.ANG'
    data_file.write_text(ANG_DATA.replace('62', '43', 1))
    assert len(EBSD(data_file).x) == 4
    with pytest.raises(ValueError, match='Unsupported EBSD file format'):
        EBSD(tmp_path/'map.unknown')


@pytest.mark.mpl_image_compare
def test_ebsd_ci():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    fig = e.plot(e.ci)
    return fig


@pytest.mark.mpl_image_compare
def test_ebsd_ci_mask():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    e.maskCI(0.1)
    fig = e.plot(e.ci)
    return fig


@pytest.mark.mpl_image_compare
def test_ebsd_ipf():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    fig = e.plotIPF()
    return fig


@pytest.mark.mpl_image_compare
def test_ebsd_ipf_1024():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    fig = e.plotIPF(1024)
    return fig


@pytest.mark.mpl_image_compare
def test_ebsd_ipf_vmask():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    e.setVMask(4)  # use only every 4th point, increases plotting speed
    e.plotIPF(1024)
    fig = e.addScaleBar()
    return fig


@pytest.mark.mpl_image_compare
def test_ebsd_ipf_crop():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    # show only a section of the image, increases plotting speed
    e.cropVMask(0, 0, 10, 10)
    fig = e.plotIPF(1024)
    return fig


@pytest.mark.mpl_image_compare
def test_ebsd_pf():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    fig = e.plotPF([1, 0, 0])
    return fig


@pytest.mark.mpl_image_compare
def test_ebsd_pf_points():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    fig = e.plotPF([1, 0, 0], points=True)
    return fig


@pytest.mark.parametrize('fileName, grid, nNeighbors', [(DATA_DIR/'EBSD.ang', 'HexGrid', 6), (None, 'SqrGrid', 4)])
def test_grid_neighbors_are_one_step_away(tmp_path, fileName, grid, nNeighbors):
    """Coordinates come from the grid; every valid neighbor is one step away."""
    if fileName is None:
        fileName = tmp_path/'square.ang'
        fileName.write_text(ANG_DATA.replace('62', '43', 1))
    ebsd = EBSD(fileName)
    assert ebsd.grid == grid
    neighbors = ebsd.neighbors()
    assert neighbors.shape == (ebsd.nPoints, nNeighbors)
    x, y = ebsd.xy()
    point, valid = np.nonzero(neighbors >= 0)
    distance = np.hypot(x[neighbors[point, valid]]-x[point], y[neighbors[point, valid]]-y[point])
    np.testing.assert_allclose(distance, ebsd.stepSizeX, rtol=1e-3)
    if grid == 'HexGrid':  # interior points have all six neighbors
        assert np.sum(np.all(neighbors >= 0, axis=1)) > 0.9*ebsd.nPoints
