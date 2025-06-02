Welcome to pythonEBSD
=====================

Electron Backscatter Diffraction (EBSD) is a microanalytical technique used in scanning electron microscopes to determine the crystallographic orientation of metals at the micrometer scale. This software package provides tools to import, analyze, and visualize the spatially resolved orientation data obtained from EBSD experiments, facilitating microstructural characterization.

Features
--------

- File formats accepted .ang | .osc | .crc | .txt
- can write .ang for FCC. Others could be added
- fast plotting interaction using virtual mask (only used for plotting)
   - increases speed in intermediate test plots
   - can be removed just before final plotting
- verified with the OIM software and mTex
- heavily tested for cubic
- separate crystal orientation and plotting of it
- some educational plotting
- examples and lots of documentation

The following modules exist
---------------------------

- :class:`ebsdlab.ebsd.EBSD`: Read EBSD data from files, plotting inverse pole figures (IPF), pole figures (PF) or normal maps
- :class:`ebsdlab.symmetry.Symmetry`: Material symmetry: cubic structure, hex, ...
- :class:`ebsdlab.orientation.Orientation`: Sum of aterial symmetry and the rotation for a material point
- :class:`ebsdlab.quaternion.Quaternion`: Mathematical description of rotations


Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   ebsd
   symmetry
   orientation
   verification
   api
