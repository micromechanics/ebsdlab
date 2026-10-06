.. _development:

Development
===========

Create an environment with Python >=3.10 and install the development requirements:

.. code-block:: console

   $ python -m pip install -r requirements-dev.txt
   $ python -m pip install .

Run the test suite from the repository root. The image tests compare plots with ``tests/baseline``, which the
documentation also shows:

.. code-block:: console

   $ pytest --mpl --mpl-baseline-path=tests/baseline

Static checks:

.. code-block:: console

   $ python -m mypy ebsdlab
   $ python -m pylint ebsdlab

Build the documentation:

.. code-block:: console

   $ make -C docs html

Open issues are tracked in `GitHub Issues <https://github.com/micromechanics/ebsdlab/issues>`_.
