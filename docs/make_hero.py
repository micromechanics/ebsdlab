"""Create docs/source/_static/hero.png (run from docs/) from EBSD.ang: raw IPF map, IPF map with unit cells, grains of
the points with CI > 0.1, IPF map after grain dilation (grain boundaries in dark grey) and the pole figure density."""
import io
import numpy as np
import scipy.ndimage as ndi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colormaps
from ebsdlab.ebsd import EBSD


def ipf(ebsd: EBSD) -> np.ndarray:
    """IPF image of the map"""
    ebsd.plotIPF('ND', show=False)
    plt.close('all')
    return ebsd.image.copy()


def addBoundaries(ebsd: EBSD, image: np.ndarray, width: int = 2) -> np.ndarray:
    """draw boundaries between grains in dark grey, width in pixels; not around points without grain"""
    grains = ebsd._image(ebsd.grainID)[0]  # pylint: disable=protected-access
    edge = np.zeros(grains.shape, dtype=bool)
    for a, b in [(grains[:, 1:], grains[:, :-1]), (grains[1:, :], grains[:-1, :])]:
        differ = (a != b) & (a > 0) & (b > 0)
        edge[:a.shape[0], :a.shape[1]] |= differ
    edge = ndi.binary_dilation(edge, np.ones((width, width))) & (grains > 0)
    image[edge] = 64
    return image


def grainColors(ebsd: EBSD) -> np.ndarray:
    """grains in the 18 non-grey tab20 colors, consecutive grains in different hues; points without grain black"""
    grains = ebsd._image(ebsd.grainID)[0]  # pylint: disable=protected-access
    palette = np.delete(colormaps['tab20'].colors, [14, 15], axis=0)  # grey is for grain boundaries
    image = (palette[grains*7 % 18]*255).astype(np.uint8)
    image[grains == 0] = 0
    return image


e = EBSD('../tests/DataFiles/EBSD.ang')
raw = ipf(e)
extent = e.imageExtent
e.maskCI(0.1)
e.removePointsOutsideMask()  # low CI points are unindexed: the band without grain is filled by dilation
e.calcGrains()
grains = addBoundaries(e, grainColors(e))
e.grainDilation()
dilated = addBoundaries(e, ipf(e))
pf = e.plotPF([1, 0, 0], show=False)
pf.axes[1].remove()  # no colorbar
buffer = io.BytesIO()
pf.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
plt.close('all')
buffer.seek(0)
pfImage = plt.imread(buffer)

fig, axes = plt.subplots(1, 5, figsize=(9, 3), dpi=120, gridspec_kw={'width_ratios': [1, 1, 1, 1, 1.9]})
for ax, image in zip(axes, [raw, raw, grains, dilated, pfImage]):
    ax.imshow(image, extent=None if image is pfImage else extent, origin='upper')
    ax.axis('off')
for x, y in [(5, 36), (16, 36), (6, 6), (9, 25)]:
    e.addUnitCellOverlay(axes[1], x, y, 3.5, 'black')
fig.tight_layout(pad=0.2, w_pad=0.5)
fig.savefig('source/_static/hero.png', facecolor='white')
