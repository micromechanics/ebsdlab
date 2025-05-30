# ebsdlab

ALPHA-VERSION: TRY AT OWN RISK

Electron Backscatter Diffraction (EBSD) is a microanalytical technique used in scanning electron microscopes to determine the crystallographic orientation of metals at the micrometer scale. This software package provides tools to import, analyze, and visualize the spatially resolved orientation data obtained from EBSD experiments, facilitating microstructural characterization.

## Features:
  - File formats accepted .ang | .osc | .crc | .txt
  - can write .ang for FCC. Others could be added
  - fast plotting interaction using virtual mask (only used for plotting)
    - increases speed in intermediate test plots
    - can be removed just before final plotting
  - verified with the OIM software and mTex
  - heavily tested for cubic
  - separate crystal orientation and plotting of it
  - some educational plotting
  - examples and lots of documentation

## Example
EBSD-Inverse Pole Figure (IPF) of polycrystalline Copper with corresponding Pole Figure
<table>
  <tr>
    <td><img src="docs/source/_static/ebsd_py_ND.png" alt="EBSD of polycrystalline Copper"></td>
    <td width="65%"><img src="docs/source/_static/ebsd_py_PF100.png" alt="Pole figure"></td>
  </tr>
</table>

## Documentation
[Documentation on github pages](https://micromechanics.github.io/ebsdlab/)

## Installation
You can install `ebsdlab` using Conda or pip.

<details>
<summary><strong>Using Conda</strong></summary>

  **Clone the repository:**

  ```console
  $ git clone https://github.com/micromechanics/ebsdlab.git ./ebsdlab
  $ cd ebsdlab
  ```

  **Create and activate the Conda environment:**

  The `environment.yml` file defines the necessary dependencies.
  ```console
  $ conda env create -f environment.yml
  ```
  After creation, activate the environment:
  ```console
  $ conda activate ebsdlab_env
  ```

  **Install the `ebsdlab` package:**
  With the Conda environment activated, install the package using pip:
  ```console
  $ python -m pip install .
  ```
</details>

<details>
<summary><strong>Using Pip</strong></summary>

  **Set up a Python environment:**
  Using a virtual environment prevents conflicts with other projects.
  ```console
  $ python -m venv venv_python_ebsd  # Create a virtual environment
  $ For Linux/macOS: source venv_python_ebsd/bin/activate
  $ For Windows: venv_python_ebsd\Scripts\activate
  ```

  **Install the `ebsdlab` package:**
  This command will install the package and dependencies:
  ```console
  $ pip install git+https://github.com/micromechanics/ebsdlab
  ```
</details>

After that, the package can be imported and used in Python codes as

```python
>>> from ebsdlab import EBSD
>>> emap = EBSD("Examples/EBSD.ang")
>>> emap.plot(e.CI)
```

## FAQ
### What features I do not envision:
  - include all crystal symmetries (materials science can mostly live with few)
  - other Euler angle definitions than Bunge; materials science does not use those

### Future features
  - improve cleaning
  - grain identification methods
  - speed up simulation
  - add different symmetries than cubic, and test

### Help wanted
 - sample files other than copper OIM files
 - feedback on tutorials
 - any feedback on functionality
 - help with cleaning and grain identification
