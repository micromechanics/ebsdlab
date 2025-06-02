.. encoding: utf-8 -*-
.. _symmetry:

Tutorials for Symmetry
======================

The **Symmetry** class is used to represent material symmetries (e.g., cubic, hex).

This example teaches the fundamentals of crystallography and shows how to determine if a vector lies within the standard stereographic triangle (SST) and retrieve its corresponding color for an Inverse Pole Figure (IPF) map.

The image shows a standard stereographic projection.

.. image:: /_static/stereographicProjection.png

Let's pick point-1 in the middle of the standard stereographic triangle: (2, 1, 3) using floats. We normalize it by dividing by its length and retrieve whether it is inside the SST and the color. We obtain that the point is inside and that the color is a very very light blue hue.

Let's pick point-2 somewhere else: (1, 2, 3). We normalize it by dividing by its length
and retrieve whether it is inside the SST and the color. We obtain that the point is outside and that the color is black, which is the default as the point is not inside the SST.

Please note, the point (2, 1, -3) is similarly inside the SST, as it is the backside projection of the sphere onto the plane.

At the end, we retrieve the standard triangle as an image.

.. jupyter-execute::

    import numpy as np
    from ebsdlab.symmetry import Symmetry
    s = Symmetry('cubic')
    point1 = np.array((2.,1.,3.))
    point1 /= np.linalg.norm(point1)
    inside, color = s.inSST( point1, color=True)
    print("The point-1 is inside:", inside)
    print("The color is:", color)
    point2 = np.array((1.,2.,3.))
    point2 /= np.linalg.norm(point2)
    inside, color = s.inSST( point2, color=True)
    print("The point-2 is inside:", inside)
    print("The color is:", color)
    s.standardTriangle()

For detailed API documentation, please refer to :py:class:`~ebsd_Symmetry.Symmetry` (material symmetry).