.. _quickstart:

Quickstart
==========

Load a map and plot the confidence index (CI), the inverse pole figure (IPF) along the normal direction (ND)
and the pole figure (PF) of the [1,0,0] direction:

- Read EBSD data from a file (.ang, .osc, .crc, .ctf, or .txt).
- Plot the confidence index.
- Plot the IPF along ND, the default, and add a scale bar.
- Plot the PF as a density.

.. jupyter-execute::

   from ebsdlab.ebsd import EBSD
   e = EBSD("../tests/DataFiles/EBSD.ang")
   print('\nPlot confidence index:')
   e.plot(e.ci)
   print('\nPlot default IPF:')
   e.plotIPF()
   print('\nSame plot as before but with scale bar:')
   e.addScaleBar()
   print('\nPlot PF as a density:')
   e.plotPF([1,0,0])

Next: the :ref:`ebsd` shows more options; :ref:`conventions` explains the directions of the plots.
