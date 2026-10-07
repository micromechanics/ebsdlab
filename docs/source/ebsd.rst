.. encoding: utf-8 -*-
.. _ebsd:

User guide
==========

Every map point has one rotation, stored in ``e.quaternions``
(``scipy.spatial.transform.Rotation``).

Data in the map
---------------

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - attribute
     - content
   * - ``e.quaternions``
     - orientation of each point
   * - ``e.x``, ``e.y``
     - coordinates from the scan grid; ``e.xy(idx)`` for selected points
   * - ``e.ci``, ``e.iq``, ``e.fit``
     - confidence index (0 bad ... 1 good), image quality, fit (EDAX)
   * - ``e.bc``
     - band contrast (Oxford)
   * - ``e.phaseID``
     - phase of each point: 0 not identified, 1, 2, ... phases; symmetry of phase k: ``e.sym[k]``
   * - ``e.grainID``
     - grain of each point: 0 no grain, 1, 2, ... grains; from ``calcGrains`` or an OIM grain file

IPF maps, pole figures, unit-cell overlays and KAM use the symmetry of each point's phase.

If a file has no readable phase symmetry, give it when opening the map (a lattice name such as ``"cubic"``,
or a ``Symmetry``; used for all phases):

.. code-block:: python

   e = EBSD("measurement.osc", symmetry="cubic")


Plotting
--------

Options that modify the plots of the :ref:`quickstart`:

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

Analysis
--------

Average orientation
~~~~~~~~~~~~~~~~~~~

Mean orientation of all points (slow; averaging over several grains is only a demonstration).

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


Grains
~~~~~~

``calcGrains`` links neighbors that are indexed, of one phase and misoriented less than ``tolerance`` (degrees);
grains are the connected regions. A grain needs at least ``minSize`` points and has to span at least ``minNRows``
rows and columns; points outside of grains have ``grainID`` 0 and are black.

``grainDilation`` then fills the points without grain, in place: a point joins the grain of most of its neighbors
and takes phase, orientation, CI and mask of the neighbor in that grain with the highest CI; this repeats until nothing
changes. It removes wild spikes and fills unindexed points; changed points are ``True`` in ``e.cleaned``. Reload
the file to undo.

.. jupyter-execute::

   import numpy as np
   from matplotlib import colormaps
   from ebsdlab.ebsd import EBSD
   cmap = colormaps["tab20"].with_extremes(bad="k")  # 20 colors, no grain in black
   e = EBSD("../tests/DataFiles/EBSD.ang")
   e.calcGrains(tolerance=5, minSize=6, minNRows=2)
   e.plot(np.where(e.grainID > 0, e.grainID % 20, np.nan), cmap=cmap, cbar=False)
   e.grainDilation()
   e.plot(np.where(e.grainID > 0, e.grainID % 20, np.nan), cmap=cmap, cbar=False)


Writing data
------------

Export to OIM
~~~~~~~~~~~~~

Remove points below a confidence index of 0.1 and write an .ang file that OIM can read again:

.. jupyter-execute::

   from ebsdlab.ebsd import EBSD
   e = EBSD("../tests/DataFiles/EBSD.ang")
   e.maskCI(0.1)
   e.removePointsOutsideMask()
   e.writeANG("ebsd.ang")

To read OIM grain files, export from OIM: Partition -> export -> grain file (type 1, saves a txt file).


Pole figures from given orientations
------------------------------------

Without a file: each ``|``-separated entry gives three Euler angles, the spread (standard deviation) around
them and the number of points.

.. jupyter-execute::

   from ebsdlab.ebsd import EBSD
   e = EBSD('void318.|125.|219.6|0.2|10')
   e.plotPF(width=10)
   e.plotPF(points=True)
