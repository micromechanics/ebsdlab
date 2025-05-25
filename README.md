Read and plot electron backscatter diffraction (EBSD) data

I look for more contributors, see below

# Features:
  - File formats accepted .ang | .osc | .crc | .txt
  - can write .ang for FCC. Others could be added
  - fast plotting interaction using virtual mask (only used for plotting)
    - increases speed in intermediate test plots
    - even though some interpolation and reading takes time
    - can be removed just before final plotting
  - verified with the OIM software and mTex
  - heavily tested for cubic
  - separate crystal orientation and plotting of it
  - some educational plotting
  - examples and lots of docummentation via Doxygen
  - doctest verifies functionality
  - introduce Kernel average misorientation


EBSD-Inverse Pole Figure (IPF) of polycrystalline Copper
![EBSD of polycrystalline Copper](docs/HTMLInputStatic/ebsd_py_ND.png)

corresponding pole-figure
![Pole figure](docs/HTMLInputStatic/ebsd_py_PF100.png)


# Documentation
[Documentation on github pages](https://steffenbrinckmann.github.io/pythonEBSD/index.html)

# What features I do not envision:
  - include all crystal symmetries (materials science can mostly live with few)
  - other Euler angle definitions than Bunge; materials science does not use those

# Future features
  - improve cleaning
  - grain identification methods
  - speed up simulation
  - add different symmetries than cubic, and test

# Help wanted
 - sample files other than copper OIM files
 - feedback on tutorials
 - any feedback on functionality
 - help with cleaning and grain identification


# Get and build docummentation
Clone data from github

```console 
$ git clone https://github.com/micromechanics/pythonEBSD.git ./pythonEBSD
```
Create conda environment and install py_ebsd package

```console
$ cd pythonEBSD
$ conda env create -f environment.yml
$ python -m pip install .
```

Build documentation and tutorials

```console
$ cd docs
$ python verifyAll.py
```
then open docs/docs/index.html in webbrowser
