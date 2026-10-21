.. encoding: utf-8 -*-
.. _comparisonOrix:

Comparison with orix and MTEX
=============================

`orix <https://orix.readthedocs.io>`_ is a Python library for crystal orientation analysis;
`MTEX for Python <https://mtex-toolbox.github.io>`_ (``pip install mtex``) is the Python version of the MATLAB
texture toolbox MTEX. This page summarizes where ebsdlab differs from them, so users can choose the right tool. The
details refer to orix 0.15.0 and MTEX for Python 1.0.0.

Advantages of ebsdlab
---------------------

- **File formats:** reads binary EDAX ``.osc`` and ``.h5`` and Oxford ``.crc``/``.cpr`` directly. orix
  reads ``.ang``, ``.ctf``, h5ebsd variants and EMsoft files, but not ``.osc``, EDAX ``.h5`` or ``.crc``.
  MTEX reads ``.osc`` and ``.crc`` too, but not the EDAX OIM ``.h5`` file of the tests (``CPTi.h5``).
- **Hexagonal grids:** handles TSL ``HexGrid`` scans natively, including the six-neighbor stencil
  for neighbor and KAM calculations. With orix, hexagonal scans are therefore usually converted to a square grid first (e.g. in EDAX
  OIM), which resamples the data. ebsdlab loads, plots and analyzes hexagonal scans directly, as does MTEX.
- **Memory:** imported Euler angles and scalar data are stored as ``float16``, which suits large
  maps; EBSD indexing accuracy (about 0.1°) does not benefit from higher input precision. ebsdlab needs the least
  memory of the three (see `Resource usage`_).
- **Lightweight:** depends only on numpy, scipy and matplotlib (see `Dependencies`_).
- **Graphical user interface:** ``ebsdlab-gui`` shows CI maps, IPF maps and pole figures, applies
  CI/crop/preview filters, places unit-cell overlays, and copies the equivalent Python code. MTEX has an import
  wizard that writes the import script; orix has no GUI.
- **Fast maps and previews:** maps are drawn directly from the scan grid, without interpolation;
  the virtual mask (``vMask``) crops the view or gives coarse previews for quick intermediate plots.
- **Built-in analysis and plotting:** KAM, grain reconstruction and dilation, unit-cell overlays on IPF maps, scale
  bars, and educational plots such as the standard triangle and unit-cell plots.
- **Verified:** results are compared against OIM and mTex for cubic materials (:ref:`verification`).
- **Simple:** small code base, Bunge Euler angles only; easy to read and to use for teaching.
- **License:** MIT; MTEX is GPL-2.0, which matters when it is built into other software.

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

Advantages of MTEX for Python
-----------------------------

- **Scope:** the full MTEX toolbox: grains and grain boundaries, denoising, parent grain reconstruction, ODFs
  from EBSD and pole figures, texture characteristics, tensors, Taylor model and optics.
- **Symmetry:** all point groups, and Miller indices, orientations and misorientations as in MATLAB MTEX.
- **Speed:** loads the ``.ang`` file and plots the IPF map faster than ebsdlab and orix (see `Resource usage`_);
  its IPF colors come from a precomputed color grid, read in a numba loop.
- **MATLAB users:** the names, syntax and documentation of MATLAB MTEX, every documentation page also in Python.


Resource usage
--------------

A one-time comparison on 2026-10-21, ebsdlab 0.0.9, orix 0.15.0 and MTEX for Python 1.0.0 (``mtex[numba,plot,io]``),
all with Python 3.14.4, numpy 2.5.4, scipy 1.18.1 and matplotlib 3.11.2, on a laptop with an AMD Ryzen 5 3500U and
14 GB memory. Each package is installed alone in a fresh virtual environment, and each measurement runs in a new
process: load the ``.ang`` file, then plot the IPF map along Z, with orix and MTEX in the TSL color key. One warm-up
run fills the numba caches of orix and MTEX; the three runs after it agree within 0.5 s. Peak memory is the peak
resident memory of the process, including the imports; loaded data is the size of all arrays of the loaded map. The
script is ``docs/benchmark.py``.

The file is ``EBSD_deformed_I_ED.ang`` of Maj et al. (Zenodo 18668585, 190 MB), 1,593,848 points. Its hexagonal
grid is made rectangular for the comparison, because orix cannot plot hexagonal grids: every second row is shifted
by half a step and the last point of the other rows is dropped. The map is then distorted, but all packages
handle the same numbers.

.. list-table::
   :header-rows: 1

   * -
     - ebsdlab
     - orix
     - MTEX
   * - load
     - 2.8 s
     - 3.1 s
     - 1.7 s
   * - IPF map
     - 2.7 s
     - 4.4 s
     - 0.9 s
   * - peak memory
     - 587 MB
     - 1235 MB
     - 916 MB
   * - memory after the imports
     - 102 MB
     - 244 MB
     - 106 MB
   * - loaded data
     - 96 MB
     - 198 MB
     - 196 MB
   * - virtual environment on disk
     - 0.30 GB
     - 0.59 GB
     - 0.54 GB

MTEX is the fastest, in loading and even more in plotting the IPF map. ebsdlab needs the least memory, about two
thirds of MTEX and half of orix, and the least disk space: the loaded map is stored as ``float16``, and the IPF
colors are computed from one pole per point, turned by the symmetry operations. On the original hexagonal file
(1,594,445 points), which orix loads but cannot plot, ebsdlab and MTEX take as long as on the rectangular one.


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
   * - **mtex[numba,plot,io] 1.0.0**
     - 6 (numpy, scipy, finufft; extras numba, matplotlib, h5py)
     - 16
