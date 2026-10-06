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


Resource usage
--------------

A one-time comparison on 2026-10-06, ebsdlab 0.0.7b1 (development version of that day) and orix 0.15.0, both with Python 3.14.4,
numpy 2.5.3, scipy 1.18.1 and matplotlib 3.11.2, on a laptop with an AMD Ryzen 5 3500U and 14 GB memory. Each
package is installed alone in a fresh virtual environment, and each measurement runs in a new process, three runs
that agree within 0.1 s: load the ``.ang`` file, then plot the IPF map along Z. Peak memory is the peak resident memory of the process,
including the imports; loaded data is the size of all arrays of the loaded map. The script is
``docs/benchmark_orix.py``.

The file is ``EBSD_deformed_I_ED.ang`` of Maj et al. (Zenodo 18668585, 190 MB), 1,593,848 points. Its hexagonal
grid is made rectangular for the comparison, because orix cannot plot hexagonal grids: every second row is shifted
by half a step and the last point of the other rows is dropped. The map is then distorted, but both packages
handle the same numbers.

.. list-table::
   :header-rows: 1

   * -
     - ebsdlab
     - orix
   * - load
     - 3.0 s
     - 3.1 s
   * - IPF map
     - 2.8 s
     - 4.2 s
   * - peak memory
     - 599 MB
     - 1235 MB
   * - memory after the imports
     - 112 MB
     - 244 MB
   * - loaded data
     - 96 MB
     - 198 MB
   * - virtual environment on disk
     - 0.30 GB
     - 0.59 GB

Both load the map equally fast; most of the time goes into parsing the text file. ebsdlab plots the IPF map faster
and needs half the memory and disk space: the loaded map is stored as ``float16``, and the IPF colors are computed
from one pole per point, turned by the symmetry operations.


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

