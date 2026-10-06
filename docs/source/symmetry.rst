.. encoding: utf-8 -*-
.. _symmetry:

Tutorial for Symmetry
=====================

The :class:`~ebsdlab.symmetry.Symmetry` class takes a lattice name: ``cubic``, ``hexagonal``, ``tetragonal``,
``orthorhombic``, ``monoclinic``, ``triclinic`` or ``trigonal`` (alias ``rhombohedral``). It holds the proper
rotations only, not the space group.

The standard stereographic triangle (SST) decides the IPF color of a direction.

.. image:: /_static/stereographicProjection.png

- point 1 = (2, 1, 3) is inside the SST: very light blue
- point 2 = (1, 2, 3) is outside: black
- (2, 1, -3) is also inside, as it is the backside of the sphere

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

``unitCell()`` returns the edges of a centered cell. The defaults only make the shape clear: cubic
``a = b = c = 1``, tetragonal and hexagonal ``c = 1.5``, orthorhombic ``a, b, c = 1, 1.25, 1.5``. Pass material
values when needed (angles in degrees; monoclinic is unique-``b``):

.. code-block:: python

    cell = Symmetry('orthorhombic').unitCell(a=2.0, b=3.0, c=4.0)

The seven crystal systems
-------------------------

.. jupyter-execute::

   import matplotlib.pyplot as plt
   from ebsdlab.symmetry import Symmetry

   lattices = (
       "cubic", "hexagonal", "tetragonal", "orthorhombic",
       "monoclinic", "triclinic", "trigonal",
   )
   cells = {lattice: Symmetry(lattice).unitCell() for lattice in lattices}
   limit = max(abs(cell[:, :]).max() for cell in cells.values()) * 1.1

   figure = plt.figure(figsize=(12, 6))
   for index, lattice in enumerate(lattices, start=1):
       axes = figure.add_subplot(2, 4, index, projection="3d")
       for edge in cells[lattice]:
           axes.plot(*edge.reshape(2, 3).T, color="C0", linewidth=2)

       axes.set_title(lattice.capitalize(), fontsize=20)
       axes.set_box_aspect((1, 1, 1))
       axes.set_xlim(-limit, limit)
       axes.set_ylim(-limit, limit)
       axes.set_zlim(-limit, limit)
       axes.set_axis_off()
       axes.view_init(elev=20, azim=35)

   figure.tight_layout()
