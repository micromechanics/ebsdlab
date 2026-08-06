.. _symmetry-gallery:

How to compare all supported crystal symmetries
================================================

Use :meth:`~ebsdlab.symmetry.Symmetry.unitCell` to obtain the centered edges
of each conventional plotting cell. The following executable example draws all
seven supported crystal systems in one Matplotlib figure. The default lengths
and angles are illustrative, so pass material-specific constants when needed.

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

``rhombohedral`` is accepted as an alias for ``trigonal``. It produces the
same rhombohedral plotting cell and proper rotational symmetry.
