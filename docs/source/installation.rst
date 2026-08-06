.. Installation

Installation
============

You can install `ebsdlab` using Conda or pip. `ebsdlab` supports Python >=3.10.

.. tabs::

   .. tab:: Using Conda

      **1. Clone the repository:**

      .. code-block:: console

         $ git clone https://github.com/micromechanics/ebsdlab.git ./ebsdlab
         $ cd ebsdlab

      **2. Create and activate the Conda environment:**

      The `environment.yml` file defines the necessary dependencies and the environment name.

      .. code-block:: console

         $ conda env create -f environment.yml

      After creation, activate the environment:

      .. code-block:: console

         $ conda activate ebsdlab

      **3. Install the `ebsdlab` package:**

      With the Conda environment activated, install the package using pip:

      .. code-block:: console

         $ python -m pip install .


   .. tab:: Using Pip

      **1. Set up a Python virtual environment (recommended):**

      Using a virtual environment prevents conflicts with other projects.

      .. code-block:: console

         $ python -m venv venv_python_ebsd  # Create a virtual environment
         #
         # Activate the environment:
         # On Linux/macOS:
         $ source venv_python_ebsd/bin/activate
         #
         # On Windows (Command Prompt):
         # venv_python_ebsd\Scripts\activate.bat
         # On Windows (PowerShell):
         # .\venv_python_ebsd\Scripts\Activate.ps1

      **2. Install the `ebsdlab` package:**

      This command will install the package and its Python dependencies directly from GitHub:

      .. code-block:: console

         $ pip install git+https://github.com/micromechanics/ebsdlab.git

After that, the package can be imported and used in Python codes as

```python
>>> from ebsdlab import EBSD
>>> emap = EBSD("tests/DataFiles/EBSD.ang")
>>> emap.plot(emap.ci)
```
