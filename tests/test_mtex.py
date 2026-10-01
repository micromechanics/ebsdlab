#!/usr/bin/python3
"""TEST EBSD class with MTEX example data; downloaded on demand since MTEX is GPL-2.0 (see DataFiles/README.md)"""
import urllib.request
import warnings
from pathlib import Path
import pytest
from ebsdlab.ebsd import EBSD


MTEX_URL = 'https://raw.githubusercontent.com/mtex-toolbox/mtex/c836b404a6729ef339857e216ff4adda143d38fb/data/EBSD/'
CACHE_DIR = Path(__file__).parent/'mtex_cache'


def mtexFile(*names: str) -> Path:
    """Download MTEX example files into the cache; skip the test with a warning if that is impossible."""
    CACHE_DIR.mkdir(exist_ok=True)
    for name in names:
        path = CACHE_DIR/name
        if path.exists():
            continue
        try:
            urllib.request.urlretrieve(MTEX_URL+name, path)
        except OSError as error:
            path.unlink(missing_ok=True)
            warnings.warn(f'MTEX example {name} cannot be downloaded: {error}')
            pytest.skip(f'MTEX example {name} cannot be downloaded')
    return CACHE_DIR/names[0]


@pytest.mark.mpl_image_compare
def test_mtex_ang_bcc_square_grid():
    e = EBSD(str(mtexFile('DC06_2uniax.ang')))
    assert e.grid == 'SqrGrid'
    return e.plotIPF(show=False)


@pytest.mark.mpl_image_compare
def test_mtex_osc_copper():
    e = EBSD(str(mtexFile('copper.osc')), symmetry='cubic')
    assert e.nPoints == 16116
    return e.plotIPF(show=False)


@pytest.mark.mpl_image_compare
def test_mtex_ang_multiphase_phase_map():
    e = EBSD(str(mtexFile('olivineopticalmap.ang')))
    return e.plot(e.phaseID, show=False)


@pytest.mark.mpl_image_compare
def test_mtex_crc_titanium_alpha_beta():
    e = EBSD(str(mtexFile('EDXLMDTi64.crc', 'EDXLMDTi64.cpr')))
    assert e.nPoints == 512*384
    return e.plotIPF(show=False)


@pytest.mark.mpl_image_compare
def test_mtex_ctf_magnesium_twins():
    e = EBSD(str(mtexFile('twins.ctf')))
    assert e.sym[1].lattice == 'hexagonal'
    return e.plotIPF(show=False)
