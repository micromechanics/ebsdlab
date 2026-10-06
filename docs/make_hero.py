"""Create docs/source/_static/hero.png (run from docs/): fcc unit cell (copper) next to the IPF map of EBSD.ang."""
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ebsdlab.ebsd import EBSD


e = EBSD('../tests/DataFiles/EBSD.ang')
e.maskCI(0.1)
e.plotIPF('ND', show=False)
image, extent = e.image, e.imageExtent
plt.close('all')

fig = plt.figure(figsize=(9, 5), dpi=150)
ax = fig.add_subplot(1, 2, 1, projection='3d')
corners = np.array(list(itertools.product((0, 1), repeat=3)), dtype=float)
for a, b in itertools.combinations(corners, 2):
    if np.isclose(np.linalg.norm(a-b), 1):
        ax.plot(*zip(a, b), color='0.25', lw=1.5)
faces = np.array([[.5, .5, 0], [.5, .5, 1], [.5, 0, .5], [.5, 1, .5], [0, .5, .5], [1, .5, .5]])
ax.scatter(*corners.T, s=420, color='#c87533', edgecolor='0.2', depthshade=False)
ax.scatter(*faces.T, s=420, color='#e8a56a', edgecolor='0.2', depthshade=False)
ax.set_box_aspect((1, 1, 1))
ax.view_init(elev=18, azim=-60)
ax.set_axis_off()
ax.set_title('Crystal: fcc unit cell', fontsize=13)

ax2 = fig.add_subplot(1, 2, 2)
ax2.imshow(image, extent=extent, origin='upper')
ax2.axis('off')
e.addUnitCellOverlay(ax2, 5, 37, 2, 'black')
e.addUnitCellOverlay(ax2, 18, 37, 2, 'black')
ax2.set_title('Map: orientation (IPF) of every point', fontsize=13)
fig.tight_layout()
fig.savefig('source/_static/hero.png', facecolor='white')
