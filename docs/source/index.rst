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

- :class:`ebsdlab.ebsd.EBSD` read/use EBSD data from files, generate random, plotting
- :class:`ebsdlab.orientation.Orientation` material symmetry + rotation of a material point
- :class:`ebsdlab.symmetry.Symmetry` material symmetry: cubic, hex, ...
- :class:`ebsdlab.quaternion.Quaternion` mathematical description of rotations
- :class:`ebsdlab.rodrigues.Rodrigues` rotation description (not used in any other file)


Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   ebsd
   verification
   orientation1
   api
