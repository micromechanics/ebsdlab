#!/usr/bin/python3
"""TEST EBSD class """
from pathlib import Path
import numpy as np
import pytest
from ebsdlab.ebsd import EBSD, SUPPORTED_SUFFIXES


DATA_DIR = Path(__file__).parent/'DataFiles'


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


@pytest.mark.mpl_image_compare
def test_ebsd_ci():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    fig = e.plot(e.CI)
    return fig


@pytest.mark.mpl_image_compare
def test_ebsd_ci_mask():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    e.maskCI(0.1)
    fig = e.plot(e.CI)
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
