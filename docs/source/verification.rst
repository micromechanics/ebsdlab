.. encoding: utf-8 -*-
.. _verification:

Verification
============

ebsdlab plots the same maps as OIM and mTex, which differ in coordinate system, projection direction and color
scheme. Here ebsdlab is compared one to one with both; maps from publications are on the next page
(:ref:`comparison-gallery`).

Same map in OIM and mTex
------------------------

Copper polycrystal
~~~~~~~~~~~~~~~~~~

.. list-table:: Comparison Table
   :header-rows: 1

   * - description
     - OIM software
     - mTex software
     - **ebsdlab**
   * - IPF ND*
     - .. image:: _static/ebsd_OIM_ND.bmp
     - .. image:: _static/ebsd_mTex_ND.png
     - .. image:: ../../tests/baseline/test_ebsd_verification_ipf_ND.png
   * - IPF RD (OIM), TD (ebsdlab)
     - .. image:: _static/ebsd_OIM_RD.bmp
     - I cannot produce
     - .. image:: ../../tests/baseline/test_ebsd_verification_ipf_TD.png

OIM's RD is the file's x-axis, which points up in the map. ebsdlab names the directions after the map
(RD = X right, TD = Y down, ND = Z into the sample; see :ref:`conventions`), so OIM's
RD is ebsdlab's -TD, which has the same IPF colors.

Issues in mTex:

- Inverse pole figure: bmp image (left) has a low color number when exporting from external window. The png figure export works better (right)

  .. image:: _static/ebsd_mTex_ND.bmp
  .. image:: _static/ebsd_mTex_ND.png

Python code to create ebsdlab results:

.. jupyter-execute::

     from ebsdlab.ebsd import EBSD
     e = EBSD("../tests/DataFiles/EBSD.ang")
     e.maskCI(0.1)
     e.plotIPF("ND")
     e.addScaleBar()
     e.plotIPF("TD")
     e.addScaleBar();


Bicrystal
~~~~~~~~~

.. list-table:: Bicrystal Comparison
   :header-rows: 1

   * - description
     - OIM software
     - mTex software
     - **ebsdlab**
   * - IPF ND*
     - .. image:: _static/bc_OIM_ND.bmp
     - .. image:: _static/bc_mTex_ND.png
     - .. image:: ../../tests/baseline/test_ebsd_bicrystal_ipf_ND.png
   * - IPF RD (OIM), TD (ebsdlab)
     - .. image:: _static/bc_OIM_RD_y.bmp
     - I cannot produce
     - .. image:: ../../tests/baseline/test_ebsd_bicrystal_ipf_TD.png
   * - PF [100]
     - .. image:: _static/bc_OIM_PF.bmp
     - .. image:: _static/bc_mTex_PF.png
     - .. image:: ../../tests/baseline/test_ebsd_bicrystal_pf.png


Python code to create the bicrystal results:

.. jupyter-execute::

     e = EBSD("../tests/DataFiles/EBSD.ang")
     e.cropVMask(ymin=35)
     e.plotIPF("ND")
     e.addSymbol(5, 37, scale=2)
     e.addSymbol(18, 37, scale=2)
     e.plotIPF("TD")
     e.addSymbol(5, 37, scale=2)
     e.addSymbol(18, 37, scale=2)
     e.plotPF([1, 0, 0], points=True)


How to run mTex
~~~~~~~~~~~~~~~

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
   % select xNorth zOutOfPlane as axis in mTex (the default gives a different pole figure)
   >> plotPDF(ebsd('Cu').orientations, Miller({1 1 1},csCopper))

Save figures from a separate window as png (the bmp color scale is broken).
