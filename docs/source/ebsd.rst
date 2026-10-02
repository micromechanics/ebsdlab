.. encoding: utf-8 -*-
.. _ebsd:

Tutorial EBSD
=============

The examples demonstrate how to read EBSD files, plot inverse pole figures, and interact with the data in the file.

The ``EBSD.quaternions`` attribute is a vector of
``scipy.spatial.transform.Rotation`` containing one rotation per map
point.


Provide phase symmetry manually
-------------------------------

When a file does
not contain readable phase-symmetry metadata, supply it while opening the map:

.. code-block:: python

   from ebsdlab.ebsd import EBSD

   e = EBSD("measurement.osc", symmetry="cubic")
   e.plotIPF()

``symmetry`` accepts a lattice name supported by ``Symmetry`` (such as
``"cubic"``), or an existing ``Symmetry``. It sets the symmetry of all phases.


Example: Read EBSD data and plot the ND IPF
-------------------------------------------

Let's looking at inverse pole figure (IPF) in normal direction (ND) and pole-figure (PF) in the [1,0,0] direction.

- Read EBSD data from a file (e.g., .ang, .osc, .crc, .ctf, or .txt).
- Plot confidence index (CI). Mask out all points with a CI less than 0.1. Initially all points are present in the mask, i.e. they are shown. By masking out points, these are removed from the mask.
- Plot inverse pole figure in normal direction
- Play with different options (e.g. 1024 pixel to see speed of plotting)
- setVMask: for fast plotting (preview: maps show every k-th row and column, pole figures every k-th point)
- Pole figure (PF) in the [1,0,0] direction. The OIM software has the top left corner as coordinate origin

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

----

The following table shows details that can improve / modify the plots.

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - Python code
     - Image
   * - .. code-block:: python

          #Plot confidence index (CI) with CI mask:
          e.maskCI(0.1)
          e.plot(e.ci)
     - .. image:: ../../tests/baseline/test_ebsd_ci_mask.png
   * - .. code-block:: python

          #Plot IPF with 1024 pixel resolution:
          e.plotIPF(1024)
     - .. image:: ../../tests/baseline/test_ebsd_ipf_1024.png
   * - .. code-block:: python

          #Plot IPF but only every 4th point
          #  increases plotting speed:
          e.setVMask(4)
          e.plotIPF(1024)
     - .. image:: ../../tests/baseline/test_ebsd_ipf_vmask.png
   * - .. code-block:: python

          #Plot section of IPF
          e.cropVMask(0,0,10,10)
          e.plotIPF(1024)
     - .. image:: ../../tests/baseline/test_ebsd_ipf_crop.png
   * - .. code-block:: python

          #Plot Pole Figure (PF) with points:
          e.plotPF([1,0,0], points=True)
     - .. image:: ../../tests/baseline/test_ebsd_pf_points.png

Data that exists and can be used for plotting in plot
-----------------------------------------------------

- OIM software:

  - e.phi1, e.phi, e.phi2 : Euler angles saved as quaternions
  - e.x, e.y : x,y coordinates, computed from the grid (e.grid, e.nRows, e.nColsOdd, e.nColsEven, e.stepSizeX, e.stepSizeY); e.xy(idx) for selected points
  - e.iq, e.ci, e.phaseID : Image Quality, confidence index (bad=0 ... good=1), phase id (0: not identified, phases from 1; symmetry of phase k: e.sym[k])
  - e.semSignal : SEM signal
  - e.fit :

- Oxford:

  - bc: band contrast


Example: Interaction with OIM Software to update grain information
------------------------------------------------------------------

.. note::
   This example is worthwhile. However, the TestB.txt file is not included in the repository. Recreate it

This example demonstrates how to process an OIM data file using this Python library and subsequently export the modified data for re-import into the OIM software.


1. How to export txt-file from OIM that can be imported into this python library:

   - Partition->export-> grain file -> use "grain file type 1" (saves a txt file)
   - Partition->export-> partition data -> save as .ang

2. Process with python: remove some masked points:

   .. jupyter-execute::

      from ebsdlab.ebsd import EBSD
      e = EBSD("../tests/DataFiles/EBSD.ang")
      e.maskCI(0.1)
      e.removePointsOutsideMask()
      e.writeANG("ebsd.ang")

3. Which can then be read in OIM again.


Example: Create artificial EBSD pattern
---------------------------------------

This example generates an artificial EBSD pattern, where parameters are separated by "|":
- The first three values specify the grain orientation as Euler angles.
- The fourth value defines the standard deviation around the mean orientation.
- The fifth value specifies the number of points to generate.

.. jupyter-execute::

   from ebsdlab.ebsd import EBSD
   e = EBSD('void318.|125.|219.6|0.2|10')
   e.plotPF(size=5)
   e.plotPF(points=True)

Example: Average orientation in file
------------------------------------

Compute the mean grain orientation for all points in the map.

   - The operation ``Orientation.average()`` is computationally intensive.
   - Averaging orientations across multiple grains is not physically meaningful, but is shown here for demonstration purposes.

.. jupyter-execute::

   import numpy as np
   from ebsdlab.orientation import Orientation
   from ebsdlab.ebsd import EBSD
   Orients = []
   e = EBSD("../tests/DataFiles/EBSD.ang")
   for i in range(e.nPoints):
       Orients.append(Orientation(quaternion=e.quaternions[i], symmetry="cubic"))
   avg = Orientation.average(Orients)
   print("Average orientation", np.round(avg.asEulers(degrees=True, standardRange=True), 0))
