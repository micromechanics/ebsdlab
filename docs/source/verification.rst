.. encoding: utf-8 -*-
.. _verification:

Compare results to that of OIM and mTex
=======================================

Trust is important for all software; comparison with other software increases trust. For EBSD analyising software this comparison is complicated by different coordinate systems, different projection directions and different color schemes. Here, we compare the results of ebsdlab with mTex and the OIM software.

Compare pythonEBSD with OIM and mTex
------------------------------------

.. list-table:: Comparison Table
   :header-rows: 1

   * - description
     - OIM software
     - mTex software
     - **ebsdlab**
   * - IPF ND*
     - .. image:: _static/ebsd_OIM_ND.bmp
          :height: 300px
     - .. image:: _static/ebsd_mTex_ND.png
          :height: 300px
     - .. image:: _static/ebsd_py_ND.png
          :height: 300px
   * - IPF RD
     - .. image:: _static/ebsd_OIM_RD.bmp
          :height: 300px
     - I cannot produce
     - .. image:: _static/ebsd_py_RD.png
          :height: 300px

Issues in mTex:

- Inverse pole figure: bmp image (left) has a low color number when exporting from external window. The png figure export works better (right)

  .. image:: _static/ebsd_mTex_ND.bmp
     :height: 300px
  .. image:: _static/ebsd_mTex_ND.png
     :height: 300px

- Pole-figure: Explicitly select x-axis as North and z-axis as outOfPlane; normal orientation has different result, although it should be the same

  .. image:: _static/ebsd_mTex_PF100Contour_xNorthzOutOfPlane.png
     :height: 200px
  .. image:: _static/ebsd_mTex_PF100Contour_org.png
     :height: 200px


Compare the three software for bicrystal
----------------------------------------

.. list-table:: Bicrystal Comparison
   :header-rows: 1

   * - description
     - OIM software
     - mTex software
     - **ebsdlab**
   * - IPF ND*
     - .. image:: _static/bc_OIM_ND.bmp
          :width: 300px
     - .. image:: _static/bc_mTex_ND.png
          :width: 300px
     - .. image:: _static/bc_py_ND.png
          :width: 300px
   * - IPF RD
     - .. image:: _static/bc_OIM_RD_y.bmp
          :width: 300px
     - I cannot produce
     - .. image:: _static/bc_py_RD.png
          :width: 300px
   * - PF [100]
     - .. image:: _static/bc_OIM_PF.bmp
          :height: 200px
     - .. image:: _static/bc_mTex_PF.png
          :height: 200px
     - .. image:: _static/bc_py_PF.png
          :height: 200px


Python code to create ebsdlab results:

.. jupyter-execute::

     from ebsdlab.ebsd import EBSD
     e = EBSD("../tests/DataFiles/EBSD.ang")
     e.cropVMask(ymin=35)
     e.plotIPF("ND")
     e.addSymbol(5,37, scale=2)
     e.addSymbol(18,37, scale=2)
     e.plotIPF("RD")
     e.addSymbol(5,37, scale=2)
     e.addSymbol(18,37, scale=2)
     e.plotIPF("TD")
     e.addSymbol(5,37, scale=2)
     e.addSymbol(18,37, scale=2)


How to run mTex
---------------

.. code-block:: matlab

   >> startup_mtex
   >> import_wizard('ebsd')
   % and select EBSD.osc
   % select plotting convention 5: x-to-right; y-to-bottom
   % select "convert Euler 2 Spatial Reference Frame"
   % save to workspace variable
   >> csCopper = ebsd('Cu').CS;
   >> plot(ebsd('Cu'),ebsd('Cu').orientations,'coordinates','on')
   >> cS = crystalShape.cube(ebsd.CS)
   >> region = [0 35 50 50];
   >> ebsdC  = ebsd(inpolygon(ebsd,region))
   >> plot(ebsdC('Cu'),ebsdC('Cu').orientations,'coordinates','on')
   >> plotPDF(ebsd('Cu').orientations, Miller({1 0 0},csCopper))
   % select xNorth zOutOfPlane as axis in mTex
   >> plotPDF(ebsd('Cu').orientations, Miller({1 1 1},csCopper))
   >> odf = calcODF(ebsd('Cu').orientations)
   >> plotPDF(odf,Miller({1 0 0},csCopper) )

If separate window: save as png, because bmp colorscale is broken
- if not separate window: save as bmp, because png crops sections off
- select xNorth zOutOfPlane as axis in mTex
- compare to original which should be the same

Example: Compare to OIM Software
--------------------------------

OIM software shows the 2D projection with the Rolling Direction (RD) upward. Note that many textbooks have the RD downward. The Normal Direction (ND) always points out of the plane; the Transverse Direction (TD) changes depending on RD.

.. jupyter-execute::

   import numpy as np
   from ebsdlab.orientation import Orientation
   o = Orientation(Eulers=np.radians([0,10,10]), symmetry="cubic")
   o.plot( )
   o.plot(plot2D='up-left')
   o.plot(poles=[1,0,0], plot2D='up-left', scale=1.5)
   o.plot(poles=[1,1,1])
   o.toScreen(equivalent=False)

Which outputs HKL and UVW as integers:
    - Euler angles: [ 0. 10. 10.]
    - HKL [ 1  5 32]
    - UVW [ 5 -1  0]

The HKL and UVW vectors are rounded to integers, hence they are approximate values. They are convenient for quick inspection but not precise.
