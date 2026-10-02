.. encoding: utf-8 -*-
.. _comparisonOrix:

Comparison with Orix
====================

`orix <https://orix.readthedocs.io>`_ is a Python library for crystal
orientation analysis. This page summarizes where ebsdlab and orix differ, so users can choose
the right tool. The orix details refer to version 0.15.0.

Advantages of ebsdlab
---------------------

- **File formats:** reads binary EDAX ``.osc`` and Oxford ``.crc``/``.cpr`` directly. orix reads
  ``.ang``, ``.ctf``, h5ebsd variants and EMsoft files, but not ``.osc`` or ``.crc``.
- **Hexagonal grids:** handles TSL ``HexGrid`` scans natively, including the six-neighbor stencil
  for neighbor and KAM calculations. With orix, hexagonal scans are therefore usually converted to a square grid first (e.g. in EDAX
  OIM), which resamples the data. ebsdlab loads, plots and analyzes hexagonal scans directly.
- **Memory:** imported Euler angles and scalar data are stored as ``float16``, which suits large
  maps; EBSD indexing accuracy (about 0.1°) does not benefit from higher input precision.
- **Lightweight:** depends only on numpy, scipy and matplotlib (see `Dependencies`_).
- **Graphical user interface:** ``ebsdlab-gui`` shows CI maps, IPF maps and pole figures, applies
  CI/crop/preview filters, places unit-cell overlays, and copies the equivalent Python code.
- **Fast maps and previews:** maps are drawn directly from the scan grid, without interpolation;
  the virtual mask (``vMask``) crops the view or gives coarse previews for quick intermediate plots.
- **Built-in analysis and plotting:** KAM, unit-cell overlays on IPF maps, scale bars, and
  educational plots such as the standard triangle and unit-cell plots.
- **Verified:** results are compared against OIM and mTex for cubic materials (:ref:`verification`).
- **Simple:** small code base, Bunge Euler angles only; easy to read and to use for teaching.

Advantages of orix
------------------

- **Symmetry:** all 32 point groups and the Laue groups, plus space groups via diffpy.structure.
  ebsdlab supports cubic, hexagonal, tetragonal and orthorhombic symmetry.
- **Phase handling:** a ``CrystalMap`` keeps a phase list with names, colors and crystal structures.
  ebsdlab also uses one symmetry per phase (via ``phaseID``), but stores only the symmetry.
- **Ecosystem:** integrates with kikuchipy (pattern indexing) and pyxem (TEM/4D-STEM). Pattern
  indexing is outside ebsdlab's scope by design.
- **Data model:** Quaternion, Rotation, Orientation and Misorientation classes, Miller indices and
  zone axes, vectorized operations on arrays of any shape, and HDF5 saving and loading.
- **Plotting:** a stereographic projection for matplotlib, IPF color keys for all Laue groups, and
  pole density plots.
- **Speed on large data:** uses numba and dask.


Dependencies
------------

.. list-table::
   :header-rows: 1

   * - package
     - direct dependencies
     - total installed
   * - **ebsdlab**
     - 3 (numpy, scipy, matplotlib)
     - 12
   * - **ebsdlab[gui]**
     - 4 (adds PySide6)
     - 16
   * - **orix 0.15.0**
     - 11
     - 41

