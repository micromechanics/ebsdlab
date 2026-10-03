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

### Group: Update documentation
- think which doc pages are out-of-date: move things around. Go away from two pages doing the same (symmetry and 7 shapes; OIM-comparison and with upublications), think about order, structure, ...
- Next: compare the examples of `docs/source/howto/comparison-gallery.rst` one by one with their sources.
  - `Ti_ZrN.ctf` vs. Kennedy et al. 2021, Fig. 10c: same field and map orientation (ZrN particle overlaps without
    flip). The paper's "IPF ND" matches ebsdlab IPF along X (RD), not Z; α-Ti also needs the 30° hex frame (X||a*):
    color distance ND 185 → RD 151 → RD with q*Rz(±30°) 94; ZrN ND 175 → RD 118. Paper maps are noise-reduced.
    Open: why X; the AZtec header says "Euler angles refer to Sample Coordinate system (CS0)". The 30° crystal
    turn is now applied by the loader.
  - `W_TKD.ctf` is published: J. Wang et al., J. Mater. Res. 37 (2022), doi:10.1557/s43578-022-00733-9 (local:
    `~/Downloads/s43578-022-00733-9.pdf`), Fig. 1: IPF RD, TD, ND of the 2000 nm wedge indent at 0.05 1/s, 17 nm
    step (original name `2000nm_005s_2TKD.ctf`). Compare all three directions; update the gallery ("no published
    image") and `tests/DataFiles/README.md`. Note: measured with Bruker e-FlashHD and plotted with MTEX, so it
    is a Bruker `.ctf`, not an Oxford reference.
  - then `Eclogite.crc` and the MTEX examples.
- `Orientation` class: still uses its own frame (ND out of plane, `plot2D='up-left'`); align it with
  `docs/source/conventions.rst`. Then update `docs/source/orientation.rst` and the end of `verification.rst`, which
  still describe the OIM layout (RD up), and point them to the conventions page.

### Group: implement new features
- Implement:grain reconstruction
- Loader for pymicro HDF5: Zenodo 12801865 (doi:10.5281/zenodo.12801865, CC-BY-4.0), CP-Ti grade 2 (hexagonal),
  `ET10_7_EBSD_post_mortem.h5` (10.4 MB; `ET10_7_EBSD_post_mortem_data_XYZ.h5`, 19.7 MB) with an OIM image of the
  same map, `ET10_7_EBSD_PM_clean_grains_OIM.tif` (0.9 MB), for comparison. The pymicro frame is unknown; it
  needs a row in `docs/source/conventions.rst`.
- Bruker `.ctf`: `docs/source/conventions.rst` has no Bruker row; every `.ctf` is treated as Oxford. Find a Bruker
  (Esprit) file with a figure made by Esprit; `W_TKD.ctf` cannot decide it alone, its paper plotted with MTEX.
- `plotPF` distribution: replace the pixel Gaussian on the stereographic image by a pole density function: von
  Mises-Fisher kernel (width in degrees) on the sphere, equal-area grid, normalized to mrd, then projected. Fixes
  rim/area distortion and gives comparable units; a step towards ODFs, which smooth in orientation space.

### Group: GUI
- Goal: extremely simple; only key parameters visible, everything else in the generated .py code.
- Flow: load once (info line) → process (min CI, crop to toolbar view, KAM on demand, cached) → plot; every change
  replots, no button. One 'Show' list: IPF (+direction), phase, IQ, CI, KAM, PF (+axis).
- Have one button at the end, render in full detail (omitting automatic vmask)
- Defaults instead of options: scale bar on, axes off, unit-cell size from map (possibly 3 settings), automatic preview vmask for large
  maps. Formats from `fileIO.LOADERS` (adds `.ctf`).
- Generated code: ample comments, full map, ends with `plt.savefig(...)`.



