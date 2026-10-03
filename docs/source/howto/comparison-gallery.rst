.. _comparison-gallery:

How to compare ebsdlab maps with their data sources
===================================================

This gallery puts the ebsdlab map of the example files next to the image published by the source of the data:
the paper that describes the dataset, or the MTEX documentation. Unlike :ref:`verification`, the comparison is
loose: the images differ in plotted direction, color key, cropping and specimen frame (ebsdlab rotates every file
into X right, Y down, Z into the sample, see :ref:`conventions`; MTEX uses its own
frame), so look for the same grains and similar colors, not for identical pixels.

The source images are downloaded when the documentation is built and are not part of ebsdlab; their licenses are
given below each example.

.. jupyter-execute::

   import urllib.request
   from io import BytesIO
   import numpy as np
   import matplotlib.pyplot as plt
   from PIL import Image
   from ebsdlab.ebsd import EBSD

   def compare(ebsd, url, title, direction='ND', turns=0):
       """Plot the ebsdlab IPF map (left) next to the image of the data source (right), if there is one.
       url: web address or local file of the source image; turns: quarter turns of the ebsdlab map, counterclockwise
       """
       ebsd.plotIPF(direction, show=False)
       plt.close()
       figure, (left, right) = plt.subplots(1, 2, figsize=(12, 5))
       left.imshow(np.rot90(ebsd.image, turns))
       left.set_title(f'ebsdlab: IPF {direction}' + (f', turned by {90*turns}°' if turns else ''))
       right.set_title(title)
       if url is not None and not url.startswith('http'):
           right.imshow(np.asarray(Image.open(url)))
       elif url is not None:
           # some servers reject Python's default user agent
           request = urllib.request.Request(url, headers={'User-Agent': 'ebsdlab-docs'})
           try:
               right.imshow(np.asarray(Image.open(BytesIO(urllib.request.urlopen(request).read()))))
           except OSError:
               right.text(0.5, 0.5, 'source image cannot be downloaded', ha='center')
       for axes in (left, right):
           axes.set_axis_off()
       figure.tight_layout()


Example files of ebsdlab
------------------------

The files are in ``tests/DataFiles``; ``tests/DataFiles/README.md`` lists their sources.

Magnesium alloy AZ31B, ``AZ31B.ang``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The source image is Fig. 4f of M. Maj, S. Musiał, M. Nowak, "Plastic work partitioning during slip- and
twinning-dominated deformation in AZ31B magnesium alloy", Metall. Mater. Trans. A,
`doi:10.1007/s11661-026-08296-8 <https://doi.org/10.1007/s11661-026-08296-8>`_,
`arXiv:2512.19548 <https://arxiv.org/abs/2512.19548>`_, CC-BY-4.0: the ⊥ ED specimen after fracture, IPF along
TD, made with EDAX OIM. The image is cropped from the arXiv version and stored in ``docs/source/_static``.
``AZ31B.ang`` keeps every 6th row and column of the original map, and the paper shows the map turned by 90°.
OIM's TD is the file's y-axis, which ebsdlab puts along -X, so ebsdlab plots RD.

The map checks the hexagonal crystal frame of EDAX files: the colors match. If the crystals were turned by 30°
about c, the blue (10-10) and green (2-1-10) grains would swap colors.

.. jupyter-execute::

   compare(EBSD('../tests/DataFiles/AZ31B.ang'), 'source/_static/AZ31B_Maj2026_Fig4f.png', 'Maj et al., Fig. 4f',
           direction='RD', turns=1)

Calcite and aragonite shell, ``Catillopecten.crc``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The source image is Fig. 6 of A. G. Checa et al., Sci. Rep. 12 (2022) 11510,
`doi:10.1038/s41598-022-15796-1 <https://doi.org/10.1038/s41598-022-15796-1>`_, CC-BY-4.0: calcite orientation
maps and pole figures of the aerials, made with Oxford Channel 5. ``Catillopecten.crc`` (Site 24) is not in
Fig. 6; the Zenodo record places the other maps in Fig. 6 and Supplementary Fig. S5.

.. jupyter-execute::

   compare(EBSD('../tests/DataFiles/Catillopecten.crc'),
           'https://media.springernature.com/full/springer-static/image/'
           'art%3A10.1038%2Fs41598-022-15796-1/MediaObjects/41598_2022_15796_Fig6_HTML.png',
           'Checa et al. 2022, Fig. 6')

Calcite aerial, ``Catillopecten_Fig6a.crc``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The map of Fig. 6a/b of the same paper (Site 8), cropped from the figure and stored in ``docs/source/_static``.
The paper colors along Z0 with its own key: 001 red, 120 green, 210 blue. The prism that Channel 5 draws in
Fig. 6b and the pole figures below the map check the frames of Oxford files: they match with the Euler angles as
stored and, for the {104} poles, the crystal turned by 30° about c.

.. jupyter-execute::

   compare(EBSD('../tests/DataFiles/Catillopecten_Fig6a.crc'), 'source/_static/Catillopecten_Checa2022_Fig6a.png',
           'Checa et al. 2022, Fig. 6a')

Eclogite, ``Eclogite.crc``
^^^^^^^^^^^^^^^^^^^^^^^^^^

The scan is a coarse point grid, which the paper shows only as pole figures: Fig. 3 of D. D. McNamara et al.,
J. Struct. Geol. (2023) 105033, `doi:10.1016/j.jsg.2023.105033 <https://doi.org/10.1016/j.jsg.2023.105033>`_,
CC-BY-4.0, with the omphacite pole figures of this scan (S6.3).

.. jupyter-execute::

   compare(EBSD('../tests/DataFiles/Eclogite.crc'),
           'https://ars.els-cdn.com/content/image/1-s2.0-S019181412300250X-gr3.jpg',
           'McNamara et al. 2023, Fig. 3')

Titanium with ZrN, ``Ti_ZrN.ctf``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The map is Fig. 10c of J. Kennedy et al., Addit. Manuf. 40 (2021) 101928,
`doi:10.1016/j.addma.2021.101928 <https://doi.org/10.1016/j.addma.2021.101928>`_. The figure is not openly
licensed, so compare with it on the publisher's page.

.. jupyter-execute::

   compare(EBSD('../tests/DataFiles/Ti_ZrN.ctf'), None, 'see Kennedy et al. 2021, Fig. 10c')

Tungsten, TKD, ``W_TKD.ctf``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

No image of this map has been published; 58 % of the points are not indexed and drawn black.

.. jupyter-execute::

   compare(EBSD('../tests/DataFiles/W_TKD.ctf'), None, 'no published image')


Example files of MTEX
---------------------

MTEX ships example files, which are not copied into ebsdlab because MTEX is licensed GPL-2.0. Download them
from the MTEX repository, as ``tests/test_mtex.py`` does, into the cache ``tests/mtex_cache``:

.. jupyter-execute::

   from pathlib import Path

   MTEX_DATA = ('https://raw.githubusercontent.com/mtex-toolbox/mtex/'
                'c836b404a6729ef339857e216ff4adda143d38fb/data/EBSD/')
   MTEX_FIGURES = 'https://mtex-toolbox.github.io/figures/'
   CACHE_DIR = Path('../tests/mtex_cache')

   def mtexFile(*names):
       """Download MTEX example files into the cache, unless they are there already."""
       CACHE_DIR.mkdir(exist_ok=True)
       for name in names:
           if not (CACHE_DIR/name).exists():
               urllib.request.urlretrieve(MTEX_DATA+name, CACHE_DIR/name)
       return str(CACHE_DIR/names[0])

The source images are from the `MTEX documentation <https://mtex-toolbox.github.io/>`_ (GPL-2.0). MTEX renumbers
them when its documentation is rebuilt, so an image may change or disappear.

Iron, square grid, ``DC06_2uniax.ang``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MTEX page `GND <https://mtex-toolbox.github.io/GND.html>`_: IPF along Y with grain boundaries.

.. jupyter-execute::

   compare(EBSD(mtexFile('DC06_2uniax.ang')), MTEX_FIGURES+'GND_01.png', 'MTEX: IPF Y', direction='TD')

Copper, ``copper.osc``
^^^^^^^^^^^^^^^^^^^^^^

MTEX page `EBSD2ODF <https://mtex-toolbox.github.io/EBSD2ODF.html>`_: IPF map.

.. jupyter-execute::

   compare(EBSD(mtexFile('copper.osc'), symmetry='cubic'), MTEX_FIGURES+'EBSD2ODF_01.png', 'MTEX: IPF')

Olivine and other minerals, ``olivineopticalmap.ang``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MTEX page `EBSDReferenceFrame <https://mtex-toolbox.github.io/EBSDReferenceFrame.html>`_: olivine IPF along Z.

.. jupyter-execute::

   compare(EBSD(mtexFile('olivineopticalmap.ang')), MTEX_FIGURES+'EBSDReferenceFrame_01.png', 'MTEX: IPF Z')

Titanium alpha and beta, ``EDXLMDTi64.crc``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MTEX page `TiBetaReconstruction <https://mtex-toolbox.github.io/TiBetaReconstruction.html>`_: alpha-Ti IPF map.

.. jupyter-execute::

   compare(EBSD(mtexFile('EDXLMDTi64.crc', 'EDXLMDTi64.cpr')), MTEX_FIGURES+'TiBetaReconstruction_01.png',
           'MTEX: IPF')

Magnesium twins, ``twins.ctf``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MTEX page `EBSDGrid <https://mtex-toolbox.github.io/EBSDGrid.html>`_: orientation map; MTEX rotates the Euler
angles by 180° about x when loading.

.. jupyter-execute::

   compare(EBSD(mtexFile('twins.ctf')), MTEX_FIGURES+'EBSDGrid_01.png', 'MTEX: orientation map')
