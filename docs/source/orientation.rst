.. default-domain:: py

Tutorial on Rotations
=====================

This tutorial introduces the fundamentals of rotations and uses the **Orientation** class as it includes the material symmetry (e.g. cubic crystal) and the rotation (e.g., rotated by 45 degrees).


Example: Specific orientations
------------------------------

This example investigates three specific orientations, the Euler angles that define them and their corresponding Inverse Pole Figure (IPF) colors for the [001] axis of the crystal. Answers are given in RGB scale:

- [100] = red   (R=1,B=0,G=0) (no rotation)
- [110] = green (R=0,B=1,G=0) (from [100] with rotation by 45 degrees around Phi)
- [111] = blue  (R=0,B=0,G=1) (from [100] with rotation of 55 and 45 degrees respectively)

Each example:

    1. creates an euler-angle-triple in radians by using 3 degree values,
    2. creates an Orientation object
    3. rounds the color for improved readability

.. jupyter-execute::

    import numpy as np
    from ebsdlab.orientation import Orientation

    angle = np.radians([0,0,0])
    o = Orientation(Eulers=angle, symmetry="cubic")
    np.round(o.IPFcolor([0,0,1]),3)

    angle = np.radians([0,45,0])
    o = Orientation(Eulers=angle, symmetry="cubic")
    np.round(o.IPFcolor([0,0,1]),3)

    angle = np.radians([0,55,45])
    o = Orientation(Eulers=angle, symmetry="cubic")
    np.round(o.IPFcolor([0,0,1]),3)

Plot unit cell
--------------

Plot unit cells and pole-figures using the orientation-class:

    - first item is one looking for
    - plot 2D projection with the coordinate system pointing up and left
    - plot also the poles and add scaling

.. jupyter-execute::

   import numpy as np
   from ebsdlab.orientation import Orientation
   angle = np.radians([0,55,45])
   o = Orientation(Eulers=angle, symmetry="cubic")
   o.toScreen()
   o.plot(plot2D='up-left')
   o.plot(poles=[1,0,0], plot2D='up-left', scale=1.5)



Example: [111] direction using vectors
--------------------------------------

The [111] direction can be challenging to define ad-hoc using Euler angles. This example demonstrates how to calculate and verify it using vectors. At the end, we plot the crystal and the poles of the [100] in the standard stereographic projection.

.. image:: /_static/ebsd_Orientation1.png

Procedure:

    1. we know the normal direction hkl = 111
    2. we pick an arbitrary direction (1,-1,0) and verify that it is perpendicular
    3. we normalize the vectors and calculate the third by the cross-product
    4. we create the rotation matrix by stacking the vectors next to eachother
    5. we verify it by plotting and obtain the Euler angles

.. jupyter-execute::

    import numpy as np
    from ebsdlab.orientation import Orientation
    hkl = np.array([1,1,1],   dtype=float)
    uvw1 = np.array([1,-1,0], dtype=float)
    print('Verify: dot product of hkl and uvw1 should be 0:', np.dot(hkl,uvw1))

    # normalize vectors and calculate uvw2
    hkl /= np.linalg.norm(hkl)
    uvw1 /= np.linalg.norm(uvw1)
    uvw2 = np.cross(hkl,uvw1)

    # create rotation matrix by stacking vectors
    rotM = np.vstack( (uvw1,uvw2,hkl) )
    print('Rotation matrix is: \n',rotM)

    # plot it and calculate Euler angles
    o = Orientation(matrix=rotM, symmetry='cubic')
    print('Euler angles are in degree: ',o.asEulers(degrees=True))
    print('The color is: ',np.round(o.IPFcolor( [0,0,1] ),3))
    o.plot()
    o.plot([1,0,0])

Example: Loop through all equivalent directions
-----------------------------------------------

This example demonstrates how to iterate through all symmetrically equivalent directions and calculate them.

  1. We create an orientation that we are interested in and a helping orientation which we use to iterate over its quaternions / directions.
  2. Obtain an equivalent crystal axis of the helper-orientation by applying a base vector (e.g., [1,0,0]).
  3. Calculate the sample direction  by transforming the sample coordinates of orientation 'o'
  4. Let's print the equivalent crystal axis and its transformed version

TODO: Understand what changed during versions, why not working

.. .. jupyter-execute::
..
.. import numpy as np
.. from ebsdlab.orientation import Orientation
.. o     = Orientation(Eulers=np.radians([0,45,0]), symmetry="cubic")
.. oHelp = Orientation(Eulers=np.array([0.,0.,0.]), symmetry="cubic")
.. for q_sym in oHelp.symmetry.symmetryQuats():
..   equivalent_crystal_axis = q_sym * (q_sym.conjugated() * np.array([1,1,0]))
..   sample_direction = o.quaternion.conjugated() * (o.quaternion * equivalent_crystal_axis)
..   print(f"Crystal Axis: {np.round(equivalent_crystal_axis, 3)}, Sample Direction: {np.round(sample_direction, 3)}")
..
.. q_sym = oHelp.symmetry.symmetryQuats()[0]
.. equivalent_crystal_axis = q_sym.conjugated() * np.array([1,0,0]) * q_sym


Example: Calculate average orientation
--------------------------------------

.. jupyter-execute::

    import numpy as np
    from ebsdlab.orientation import Orientation
    a = Orientation(Eulers=np.radians([0,45,0]), symmetry='cubic')
    b = Orientation(Eulers=np.radians([0,0,0]),  symmetry='cubic')
    c = Orientation(Eulers=np.radians([0,15,0]), symmetry='cubic')
    avg = Orientation.average([a,b,c])
    print("Rotation angles",avg.asEulers(degrees=True))
