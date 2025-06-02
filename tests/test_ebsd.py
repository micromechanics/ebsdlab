#!/usr/bin/python3
"""TEST EBSD class """
from pathlib import Path
import pytest
from ebsdlab.ebsd import EBSD

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
    fig =  e.plotIPF(1024)
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
    e.cropVMask(0,0,10,10)  # show only a section of the image, increases plotting speed
    fig = e.plotIPF(1024)
    return fig

@pytest.mark.mpl_image_compare
def test_ebsd_pf():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    fig = e.plotPF([1,0,0])
    return fig

@pytest.mark.mpl_image_compare
def test_ebsd_pf_points():
    dataDir = Path(__file__).parent/'DataFiles'
    e = EBSD(str(dataDir/'EBSD.ang'))
    fig = e.plotPF([1,0,0], points=True)
    return fig
