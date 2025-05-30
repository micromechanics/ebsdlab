Example: Compare with OIM software and verify data (TODO)
---------------------------------------------------------

.. code-block:: python

   import numpy as np
   from ebsdlab.ebsd import EBSD
   from ebsdlab.orientation import Orientation
   e = EBSD("Examples/EBSD.ang", doctest=True)
   # Load .ang file:  Examples/EBSD.ang
   #    Read file with step size: 0.2 0.2
   #    Optimal image pixel size: 103
   #    Number of points: 23909
   e.maskCI(0.001)
   e.plotIPF('ND')  # doctest: +SKIP
   e.cropVMask(xmin=18, ymin=12, ymax=17)
   e.plotIPF('ND')  # doctest: +SKIP

Inspect the original data at x,y = 20.1,14.38 um (line 8785 from EBSD.ang):

   5.55763   2.18448   3.83349     20.10000     14.37602 3653.496  0.482  0      1  1.104

.. code-block:: python

   print(np.degrees([5.55763, 2.18448, 3.83349]))
   # [318.4287431  125.16148443 219.64279781]
   e.cropVMask(xmin=20, xmax=20.2, ymin=14.3, ymax=14.4)
   e.y[e.vMask]  # verify y: correct if rounding accounted for
   # array([14.37602])
   angle = e.quaternions[e.vMask].asEulers().flatten()
   print(np.round(np.degrees(angle)))  # convert to only positive values
   # [ -42.  125. -140.]
   print(np.round(np.degrees(angle)+np.array([360,0,360])))
   # [318. 125. 220.]

Plot correct unit cells and pole-figures using the orientation-class:

.. code-block:: python

   o = Orientation(Eulers=angle, symmetry="cubic")
   o.doctest = True
   o.toScreen()  # first item is one looking for #doctest: +SKIP
   o.plot(plot2D='up-left')
   doctestImage("ebsd_10")
   o.plot(poles=[1,0,0], plot2D='up-left', scale=1.5)
   doctestImage("ebsd_11")
