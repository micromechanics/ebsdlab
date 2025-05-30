.. Installation

Installation
============

You can install `pythonEBSD` using Conda or pip.

.. tabs::

   .. tab:: Using Conda

      **1. Clone the repository:**

      .. code-block:: console

         $ git clone https://github.com/micromechanics/pythonEBSD.git ./pythonEBSD
         $ cd pythonEBSD

      **2. Create and activate the Conda environment:**

      The `environment.yml` file defines the necessary dependencies and the environment name (e.g., `pythonEBSD_env`).

      .. code-block:: console

         $ conda env create -f environment.yml

      After creation, activate the environment:

      .. code-block:: console

         $ conda activate pythonEBSD_env

      **3. Install the `pythonEBSD` package:**

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

      **2. Install the `pythonEBSD` package:**

      This command will install the package and its Python dependencies directly from GitHub:

      .. code-block:: console

         $ pip install git+https://github.com/micromechanics/pythonEBSD.git

After that, the package can be imported and used in Python codes as

```python
>>> import pythonEBSD
>>> emap = pythonEBSD.EBSD("Examples/EBSD.ang")
>>> emap.plot(e.CI)
```
