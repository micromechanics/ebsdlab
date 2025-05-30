.. encoding: utf-8 -*-
.. _verification:

Verification
============


Compare pythonEBSD with OIM and mTex (TODO)
-------------------------------------------

.. jupyter-execute::

   from ebsdlab.ebsd import EBSD
   e = EBSD("../tests/DataFiles/EBSD.ang")
   e.maskCI(0.001)
   e.plotPF(size=1)
   e.plotPF()
   e.plotPF(proj2D='down-right')

.. list-table:: Comparison Table
   :header-rows: 1

   * - description
     - OIM software
     - mTex software
     - this python code
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

Issues from mTex:

- .bmp (left): low color number when exporting from external window. .png works (right)

  .. image:: _static/ebsd_mTex_ND.bmp
     :height: 300px
  .. image:: _static/ebsd_mTex_ND.png
     :height: 300px

- Explicitly select x-axis as North and z-axis as outOfPlane; normal orientation has different result, although it should be the same

  .. image:: _static/ebsd_mTex_PF100Contour_xNorthzOutOfPlane.png
     :height: 200px
  .. image:: _static/ebsd_mTex_PF100Contour_org.png
     :height: 200px


Compare the three software for bicrystal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Bicrystal Comparison
   :header-rows: 1

   * - description
     - OIM software
     - mTex software
     - this python code
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

Python code:

.. code-block:: python

   from ebsdlab.ebsd import EBSD
   e = EBSD("../tests/DataFiles/EBSD.ang")
   e.cropVMask(ymin=35)
   e.plotIPF("ND", fileName="doctest.png")
   e.addSymbol(5,37, scale=2)
   e.addSymbol(18,37, scale=2, fileName="doctest.png")
   e.plotIPF("RD", fileName="pythonRD.png")
   e.addSymbol(5,37, scale=2)
   e.addSymbol(18,37, scale=2, fileName="pythonRD.png")
   e.plotIPF("TD", fileName="pythonTD.png")
   e.addSymbol(5,37, scale=2)
   e.addSymbol(18,37, scale=2, fileName="pythonTD.png")
   e.plotPF(fileName="pythonPF.png")

.. note::
   The last image is not colored. This is not implemented yet. TODO

How to run mTex
---------------

.. code-block:: matlab

   >> startup_mtex
   >> import_wizard('ebsd')
   % and select EBSD.osc
   % select plotting convention 5: x-to-right; y-to-bottom
   % select "convert Euler 2 Spacial Referecence Frame"
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

Data that exists and can be used for plotting in plot:
------------------------------------------------------

- OIM software:
  - e.phi1, e.PHI, e.phi2 : Euler angles saved as quaternions
  - e.x, e.y : x,y coordinates
  - e.IQ, e.CI, e.phaseID : Image Quality, confidence index (bad=0 ... good=1), phase id
  - e.SEMsignal : SEM signal
  - e.fit :
- Oxford:
  - bc: band contrast

Hints for developers
--------------------

- run ``./verifyAll.py`` after all changes to verify the code and create the html-documentation
- git commands:

  .. code-block:: bash

     git add -A
     git gui
     git commit -m "solved symbolic link issue"
     git push -u origin master


