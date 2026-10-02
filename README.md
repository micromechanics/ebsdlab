# ebsdlab

Electron Backscatter Diffraction (EBSD) is a microanalytical technique used in scanning electron microscopes to determine the crystallographic orientation at the micrometer scale. This software package provides tools to import, analyze, and visualize the data.

## Features:
  - File formats accepted .ang | .osc | .crc | .ctf | .txt
  - fast plotting: maps are drawn directly from the scan grid
    - virtual mask (only used for plotting)
      - increases speed in intermediate test plots
      - can be removed just before final plotting
  - verified with the OIM software and mTex
  - separate crystal orientation and plotting of it
  - some educational plotting
  - examples and documentation
  - minimal requirements on libraries


## Example
EBSD-Inverse Pole Figure (IPF) of polycrystalline Copper with corresponding Pole Figure
<table>
  <tr>
    <td><img src="docs/source/_static/ebsd_py_ND.png" alt="EBSD of polycrystalline Copper"></td>
    <td width="65%"><img src="docs/source/_static/ebsd_py_PF100.png" alt="Pole figure"></td>
  </tr>
</table>


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
  $ conda activate ebsdlab
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

After that, the package can be used as

```python
>>> from ebsdlab import EBSD
>>> emap = EBSD("tests/DataFiles/EBSD.ang")
>>> emap.plot(emap.ci)
```

### Graphical user interface

Install the optional GUI dependency and start the application:

```console
$ pip install 'ebsdlab[gui]'
$ ebsdlab-gui
```

## Documentation
[Documentation on github pages](https://micromechanics.github.io/ebsdlab/)

## FAQ
### What features I do not envision:
  - include all crystal symmetries (materials science can mostly live with few)
  - other Euler angle definitions than Bunge; materials science does not use those

### Future features
  - improve cleaning
  - grain identification methods
  - speed up simulation
  - test non-cubic symmetries further: example data covers hexagonal, trigonal, orthorhombic, monoclinic and
    triclinic phases; only hexagonal has image tests

### Help wanted
 - sample files with tetragonal phases
 - feedback on tutorials
 - any feedback on functionality
 - help with cleaning and grain identification


## Design scope
`ebsdlab` analyzes already indexed EBSD orientation data for known phases. It is
not intended to index raw Kikuchi patterns, identify phases from diffraction
patterns, or simulate EBSD patterns.

Accordingly, it models the rotational crystal symmetry required for orientation
analysis, not complete atomistic crystal structures. Atomic basis positions,
lattice centering, structure factors, and translational space-group operations
such as glide planes and screw axes are out of scope. These details are  not important
for orientation maps, IPF colors, misorientation, KAM, or pole-figure analysis of a known phase.

## Phases
Every point has a `phaseID` (`uint8`):
- `0`: not identified; `emap.sym[0]` is an empty `Symmetry()`
- `1, 2, ...`: phases; `emap.sym[k]` is the symmetry of phase `k`

IPF maps, pole figures, unit-cell overlays and KAM use the symmetry of each point's phase; KAM ignores
neighbors of another phase.

## Requirements
`ebsdlab` supports Python >=3.10.


## Development
Use the local `.venv/` when it has been prepared for this repository, or create an environment with Python >=3.10 and install the development requirements.

```console
$ python -m pip install -r requirements-dev.txt
$ python -m pip install .
```

Run the test suite from the repository root:

```console
$ pytest --mpl --mpl-baseline-path=tests/baseline
```

Static checks are configured for:

```console
$ python -m mypy ebsdlab
$ python -m pylint ebsdlab
```

## Issues
Open issues are tracked in [GitHub Issues](https://github.com/micromechanics/ebsdlab/issues). Local issue notes may also be documented in this repository when they need to stay alongside the code.

### Group: check orientations
- ensure that all doc image-files from ebsdlab come from builds;
- add more checks to compare other file outputs with those from references/Orix...
- Euler angles of all formats are used as stored. Hex/trigonal crystal frame: `.ctf`, `.crc` use X||a*, ebsdlab and
  `.ang`, `.osc` X||a (30° about c, fix: q*Rz(-30°)). Specimen frame: MTEX rotates by 180° about [1-10] (`.ang`
  setting 2) or Z (`.ctf`, `.crc`) into the map frame; ebsdlab does not, so PFs and RD/TD IPFs differ by vendor.

### Group: implement new features
- Implement:grain reconstruction
- `plotPF` distribution: replace the pixel Gaussian on the stereographic image by a pole density function: von
  Mises-Fisher kernel (width in degrees) on the sphere, equal-area grid, normalized to mrd, then projected. Fixes
  rim/area distortion and gives comparable units; a step towards ODFs, which smooth in orientation space.
- `plotPF`: `xs, ys = [], []` is inside the phase loop, so a multi-phase pole figure shows only the last phase;
  move it above the loop.

### Group: GUI
- Goal: extremely simple; only key parameters visible, everything else in the generated .py code.
- Flow: load once (info line) → process (min CI, crop to toolbar view, KAM on demand, cached) → plot; every change
  replots, no button. One 'Show' list: IPF (+direction), phase, IQ, CI, KAM, PF (+axis).
- Have one button at the end, render in full detail (omitting automatic vmask)
- Defaults instead of options: scale bar on, axes off, unit-cell size from map (possibly 3 settings), automatic preview vmask for large
  maps. Formats from `fileIO.LOADERS` (adds `.ctf`).
- Generated code: ample comments, full map, ends with `plt.savefig(...)`.



