# from ebsdlab.orientation import Orientation
# from ebsdlab.ebsd import EBSD
# import numpy as np
# Example: Verify numbers
# -----------------------

# Verify the data by using ebsdlab and OIM software.

# .. code-block:: python

# e = EBSD('Examples/EBSD.ang')
# e.cropVMask(xmin=18, ymin=12, ymax=17)
# e.plotIPF('ND')

# Inspect the original data at x, y = 20.1, 14.38 um(line 8785 from EBSD.ang):

#     5.55763   2.18448   3.83349     20.10000     14.37602 3653.496  0.482  0      1  1.104

# .. code-block:: python

# print(np.degrees([5.55763, 2.18448, 3.83349]))
# # [318.4287431  125.16148443 219.64279781]
# e.cropVMask(xmin=20, xmax=20.2, ymin=14.3, ymax=14.4)
# e.y[e.vMask]  # verify y: correct if rounding accounted for
# # array([14.37602])
# angle = e.quaternions[e.vMask].asEulers().flatten()
# print(np.round(np.degrees(angle)))  # convert to only positive values
# # [ -42.  125. -140.]
# print(np.round(np.degrees(angle)+np.array([360, 0, 360])))
# # [318. 125. 220.]
