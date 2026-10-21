.. encoding: utf-8 -*-
.. _comparisonOrix:
.. _comparisonLibraries:

Comparison with orix and MTEX
=============================

`orix <https://orix.readthedocs.io>`_ is a Python library for crystal orientation analysis;
`MTEX for Python <https://mtex-toolbox.github.io>`_ (``pip install mtex``) is the Python version of the MATLAB
texture toolbox MTEX. This page compares them with ebsdlab, so users can choose the right tool. The details refer
to orix 0.15.0 and MTEX for Python 1.0.0.

Features
--------

.. list-table::
   :header-rows: 1
   :widths: 16 28 28 28

   * -
     - ebsdlab
     - orix
     - MTEX for Python
   * - file formats
     - ``.ang``, ``.ctf``, ``.osc``, ``.crc``/``.cpr``, EDAX ``.h5``, ``.txt``
     - ``.ang``, ``.ctf``, Bruker and EMsoft h5ebsd, orix HDF5
     - ``.ang``, ``.ctf``, ``.osc``, ``.crc``/``.cpr``, HDF5 of Bruker, EDAX, Oxford, Thermo Fisher and EMSphInx;
       not EDAX ``.h5`` written by OIM 7.3 (keys start with a space), e.g. ``CPTi.h5`` of the tests
   * - hexagonal scan grids
     - load, maps and neighbor analysis
     - load only, no maps
     - load, maps and neighbor analysis
   * - crystal symmetry
     - Laue groups of the seven crystal systems
     - all 32 point groups, space groups via diffpy.structure
     - all 32 point groups
   * - maps
     - IPF, CI, IQ and other scalar maps, unit-cell overlays, scale bars
     - IPF and property maps
     - IPF and property maps, many color keys
   * - pole figures and texture
     - pole figures with the pole density in multiples of random (mrd)
     - pole figures and pole density functions
     - pole figures, ODFs from EBSD and from pole-figure data, texture characteristics
   * - KAM and grains
     - KAM, grain reconstruction, grain dilation
     - –
     - KAM, grains and grain boundaries, denoising, parent grain reconstruction
   * - graphical user interface
     - ``ebsdlab-gui``: maps, pole figures, filters, unit-cell overlays, copies the Python code
     - –
     - import wizard that writes the import script
   * - speed and memory
     - see `Resource usage`_
     - see `Resource usage`_
     - see `Resource usage`_
   * - installed packages
     - 12 (see `Dependencies`_)
     - 41
     - 4, or 16 with the numba, plot and io extras

Choose ebsdlab for large maps on little memory, hexagonal scans and quick looks: imported Euler angles and scalar
data are stored as ``float16`` (EBSD indexing accuracy, about 0.1°, does not benefit from higher precision), maps are
drawn directly from the scan grid without interpolation, and the virtual mask (``vMask``) crops the view or gives
coarse previews. Its small code base with Bunge Euler angles only is easy to read and suits teaching, and its
results are compared against OIM and MTEX (:ref:`verification`).

Choose orix for its data model: Quaternion, Rotation, Orientation and Misorientation classes, Miller indices and zone
axes, vectorized operations on arrays of any shape, and all point groups and space groups. A ``CrystalMap`` keeps a
phase list with names, colors and crystal structures; some computations run lazily with dask.

Choose MTEX for Python for the widest analysis: grains and boundaries, ODFs and texture characteristics, parent grain
reconstruction, all with the names, syntax and documentation of MATLAB MTEX. Once numba has compiled its functions
it is the fastest of the three; the first run after an installation is the slowest.


Resource usage
--------------

A one-time comparison on 2026-10-21, ebsdlab 0.0.9, orix 0.15.0 and MTEX for Python 1.0.0 (``mtex[numba,plot,io]``),
all with Python 3.14.4, numpy 2.5.4, scipy 1.18.1 and matplotlib 3.11.2, on a laptop with an AMD Ryzen 5 3500U and
14 GB memory. Each package is installed alone in a fresh virtual environment, and each measurement runs in a new
process: load the ``.ang`` file, then plot the IPF map along Z, with orix and MTEX in the TSL color key. orix and
MTEX compile numba functions on first use and cache them; their times read "first run / later runs", the first run
starting from an empty numba cache as after an installation. Three repetitions agree within 0.2 s. Peak memory is
the peak resident memory of the process, including the imports; loaded data is the size of all arrays of the loaded
map. The script is ``docs/benchmark.py``.

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
     - 2.7 s
     - 3.0 / 3.0 s
     - 5.2 / 1.3 s
   * - IPF map
     - 1.0 s
     - 4.1 / 4.0 s
     - 5.1 / 0.7 s
   * - peak memory
     - 569 MB
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

Once compiled, MTEX is the fastest: its numba parser reads the text file in half the time of ``np.loadtxt``, which
ebsdlab and orix use. On the first run after an installation, the compilation makes MTEX the slowest, 10 s for loading
and plotting together. orix compiles little. ebsdlab needs no compilation, the least memory (about 60 % of MTEX and
half of orix) and the least disk space: the loaded map is stored as ``float16``, and each pole is folded directly into
the standard triangle and colored once. On the original hexagonal file (1,594,445 points), which orix loads but cannot
plot, ebsdlab and MTEX take as long as on the rectangular one.


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
   * - **mtex 1.0.0**
     - 3 (numpy, scipy, finufft)
     - 4
   * - **mtex[numba,plot,io] 1.0.0**
     - 6 (numpy, scipy, finufft; extras numba, matplotlib, h5py)
     - 16
