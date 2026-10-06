.. _conventions:

Conventions
===========

EBSD vendors and programs use different coordinate systems for the same data. ebsdlab uses one set of
conventions for everything it computes and plots; every loader's result is rotated into them. This page collects
them in one place.

Sample frame: X, Y, Z
---------------------

The sample frame follows the map as it is shown:

- **X** points right, **Y** points down, **Z** points into the sample (right-handed: X × Y = Z)
- RD = X, TD = Y, ND = Z
- pixel (0, 0) is at the top left
- the sample frame is also called the specimen frame or sample coordinate system; it is written in upper case
- pole figures and the unit cells drawn on maps are oriented the same way: X right, Y down, seen from above the
  sample (upper hemisphere, toward -Z)

Crystal frame: x, y, z
----------------------

The crystal frame is fixed to the lattice and is written in lower case.

- cubic: x = [100], y = [010], z = [001]; all vendors agree
- hexagonal and trigonal: the lattice axes a1, a2 (120° apart) and c are not perpendicular, so the orthonormal
  axes must be chosen, and vendors choose differently. ebsdlab: x = a1 = [2-1-10], z = c, y = z × x (90° from a1,
  30° before a2)

Orientations
------------

An orientation q rotates the crystal frame into the sample frame: ``q.apply(v_crystal) = v_sample``. Euler angles
are Bunge angles (ZXZ, intrinsic, active), as in ``scipy.spatial.transform.Rotation.from_euler('ZXZ', ...)``.

Files and vendors
-----------------

After a loader has read a file, ``fileIO.rotateToConventions`` rotates the orientations into the conventions above:
q = SAMPLE · q_file · CRYSTAL. ``writeANG`` rotates back into the EDAX frame. The rotations are in
``fileIO.FRAMES``.

.. list-table::
   :header-rows: 1

   * - vendor
     - files
     - sample frame (SAMPLE)
     - crystal frame, hexagonal/trigonal (CRYSTAL)
   * - EDAX (TSL OIM)
     - ``.ang`` ``.osc`` ``.h5`` ``.txt``
     - 180° about [1-10]: file x → -Y, y → -X, z → -Z
     - none, file x = a1
   * - Oxford (AZtec, Channel 5)
     - ``.ctf`` ``.crc``
     - none: file x = X, y = Y, z = Z
     - 30° about c, file x = a1*

``.txt`` is treated as the EDAX text export, saved like ``.ang``.

How the conventions were checked
--------------------------------

- **EDAX, sample frame:** the cube symbols that OIM draws on the bicrystal of ``EBSD.ang`` match ebsdlab's only
  for 180° about [1-10] (:ref:`verification`). As a consequence, OIM's RD is ebsdlab's -TD.
- **EDAX, crystal frame:** the IPF colors of ``AZ31B.ang`` match OIM's in Maj et al., Fig. 4f; turning the
  crystals by 30° about c swaps blue and green (:ref:`comparison-gallery`).
- **Oxford, sample frame:** ``Catillopecten.crc`` (Checa et al. 2022, Fig. S6b, Channel 5), ``Eclogite_Fig5.crc``
  (McNamara et al. 2023, Fig. 5, Channel 5) and ``Ti_ZrN.ctf`` (Kennedy et al. 2021, Fig. 10c, AZtec) match their
  published maps without rotation or flip (:ref:`comparison-gallery`).
- **Oxford, crystal frame:** the α-Ti colors of ``Ti_ZrN.ctf`` match AZtec's only with the crystal turned by 30°
  about c: mean color distance 85, without the turn 151. The {104} poles that Channel 5 draws for
  ``Catillopecten_Fig6a.crc`` in Checa et al. 2022, Fig. 6, agree.
