.. _scope:

Design scope
============

ebsdlab analyzes already indexed EBSD orientation data for known phases. It is not intended to index raw
Kikuchi patterns, identify phases from diffraction patterns, or simulate EBSD patterns.

Accordingly, it models the rotational crystal symmetry required for orientation analysis, not complete atomistic
crystal structures. Atomic basis positions, lattice centering, structure factors, and translational space-group
operations such as glide planes and screw axes are out of scope. These details are not important for
orientation maps, IPF colors, misorientation, KAM, or pole-figure analysis of a known phase.

Also not planned:

- all crystal symmetries: materials science can mostly live with few (see :ref:`symmetry`)
- Euler angle definitions other than Bunge; materials science does not use those
