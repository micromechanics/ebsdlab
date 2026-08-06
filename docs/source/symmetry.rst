.. encoding: utf-8 -*-
.. _symmetry:

Tutorial for Symmetry
=====================

The :class:`~ebsdlab.symmetry.Symmetry` class uses lattice names:
``cubic``, ``hexagonal``, ``tetragonal``, and ``orthorhombic``. Internally,
SciPy represents their proper rotational symmetries by the point-group codes
``O``, ``D6``, ``D4``, and ``D2``. These SciPy groups are rotation groups, not crystallographic
space groups.

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

Unit cells
----------

``unitCell()`` returns the edges of a centered cell for every supported lattice.
Its defaults are chosen to make the lattice shape clear in plots;
they are not material-specific lattice parameters, as those differ from one material to the next. Cubic cells use
``a = b = c = 1``; tetragonal and hexagonal cells use ``a = b = 1`` and
``c = 1.5``; orthorhombic cells use ``a = 1``, ``b = 1.25``, and ``c = 1.5``.
Positive lattice constants can be supplied when required:

.. code-block:: python

    cell = Symmetry('orthorhombic').unitCell(a=2.0, b=3.0, c=4.0)

For detailed API documentation, refer to
:class:`~ebsdlab.symmetry.Symmetry`.
