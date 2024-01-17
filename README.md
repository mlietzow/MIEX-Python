# MIEX for Python

[![arXiv](https://img.shields.io/badge/arXiv-astro--ph%2F0406118-b31b1b)](https://arxiv.org/abs/astro-ph/0406118)
[![ascl](https://img.shields.io/badge/ascl-1810.019-262255)](https://ascl.net/1810.019)
[![bibcode](https://img.shields.io/badge/bibcode-2004CoPhC.162..113W-1c459b)](https://ui.adsabs.harvard.edu/abs/2004CoPhC.162..113W)
[![doi](https://img.shields.io/badge/doi-10.1016%2Fj.cpc.2004.06.070-fab70c)](https://doi.org/10.1016/j.cpc.2004.06.070)

is a Mie scattering code for large grains written in Python and based on **miex** [(Wolf & Voshchinnikov 2004)](https://ui.adsabs.harvard.edu/abs/2004CoPhC.162..113W).

The following quantities for

1. single grain sizes / chemical components and
2. mixtures of chemically different grains with a size distributions

can be calulated:

- Scattering matrix elements $S_{11}$, $S_{12}$, $S_{33}$, and $S_{34}$,
- Extinction effiency factor ($Q_\mathrm{ext}$) and Extinction cross section ($C_\mathrm{ext}$),
- Scattering effiency factor ($Q_\mathrm{sca}$) and Scattering cross section ($C_\mathrm{sca}$),
- Absorption effiency factor ($Q_\mathrm{abs}$) and Absorption cross section ($C_\mathrm{abs}$),
- Backscattering effiency factor ($Q_\mathrm{bk}$) and Backscattering cross section ($C_\mathrm{bk}$),
- Radiation pressure effiency factor ($Q_\mathrm{pr}$),
- Albedo,
- Scattering assymetry factor ($g$).


## Requirements

To run miex for Python, at least **Python 3.6** and the following packages are required:
 - numba
 - numpy

To check the current version or to install the packages, run:

```bash
python3 --version
pip3 install -r requirements.txt
```

## Run miex

[run_miex.py](run_miex.py) executes miex with the parameters of the input file:

```bash
python3 run_miex.py example1.input
python3 run_miex.py example2.input
```

The results are stored in the `results` directory.

The input file is organized as follows:

| Parameter                                                                      | Type  |
| ------------------------------------------------------------------------------ | ----- |
| Real refractive index of the surrounding medium                                | float |
| Number of wavelengths                                                          | int   |
| Number of components                                                           | int   |
| Name of the dust data files (lambda/n/k data)                                  |       |
| &ensp; 1. component                                                            | str   |
| &ensp; ...                                                                     | str   |
| &ensp; n. component                                                            | str   |
| Relative abundances of the different components [%] <sup>1</sup>               |       |
| &ensp; 1. component                                                            | float |
| &ensp; ...                                                                     | float |
| &ensp; n. component                                                            | float |
| Single grain size (1) or grain size distribution (2)                           | int   |
| Grain size [micron] / Minimum grain size [micron]                              | float |
| Maximum grain size [micron] <sup>2</sup>                                       | float |
| Size distribution exponent <sup>2</sup>                                        | float |
| Number of size bins <sup>2</sup>                                               | int   |
| Calculate scattering matrix elements (0: no / 1: yes)                          | int   |
| Number of scattering angles in the interval [0, 180]; odd number! <sup>3</sup> | int   |
| Project name for the output files                                              | str   |
| Save results in separate files (0: no / 1: yes)                                | int   |

<sup>1</sup> only for multiple compositions, omit for single composition;

<sup>2</sup> only for a grain size distribution, omit for single grain size;

<sup>3</sup> only if scattering matrix elements are calculated, omit if not.


### Run miex via streamlit

[miex_app.py](miex_app.py) can be used as a simple Streamlit web app. Install streamlit

```bash
pip3 install streamlit
```

and start the web app:

```bash
streamlit run miex_app.py
```


## Import miex

[miex.py](miex.py) can be imported and used in any python script (see [example script](miex_example.py) or [jupyter notebook](miex_notebook.ipynb)):

```bash
python3 miex_example.py
```

To calculate the efficiency factors and scattering amplitude functions (optionally), use e.g.,

```python
miex.shexqnn2(x=1.0, m=complex(1.5, 0.0), nang=91, doSA=True)
```

| Variable     | Input Parameter                                           | Type           |
| ------------ | --------------------------------------------------------- | -------------- |
| `x`          | Size parameter                                            | float          |
| `m`          | Complex refractive index                                  | complex        |
| `nang=1`     | Half number of scattering angles in the intervall [0, 90] | int, optional  |
| `doSA=False` | Calculate scattering amplitude functions                  | bool, optional |

| Variable      | Output Parameter              | Type                | Index |
| ------------- | ------------------------------| ------------------- | ----- |
| $Q_{\rm ext}$ | Extinction efficiency         | float               | 0     |
| $Q_{\rm abs}$ | Absorption efficiency         | float               | 1     |
| $Q_{\rm sca}$ | Scattering efficiency         | float               | 2     |
| $Q_{\rm bk}$  | Backscattering effiency       | float               | 3     |
| $Q_{\rm pr}$  | Radiation pressure effiency   | float               | 4     |
| $A$           | Single scattering albedo      | float               | 5     |
| $g_{\rm sca}$ | Scattering assymetry factor   | float               | 6     |
| $S_{1}$       | Scattering amplitude function | complex, array-like | 7     |
| $S_{2}$       | Scattering amplitude function | complex, array-like | 8     |

The scattering amplitude functions are an array with size `2*nang-1`.

To calculate the scattering matrix elements, use

```python
miex.scattering_matrix_elements(S1, S2)
```

where $S_1$ and $S_2$​ are the scattering amplitude functions.
The function returns the scattering matrix elements $S_{11}$​, $S_{12}$​, $S_{33}$​, $S_{34}$​ (in this order).


## Test miex

[test_miex.py](test_miex.py) includes some test routines. The results are compared with results by [Bohren & Huffman (1998)](https://doi.org/10.1002/9783527618156) and by [Wiscombe (1979)](https://doi.org/10.5065/D6ZP4414):

```bash
python3 test_miex.py
```


## Project structure

    .
    ├── ri-data                                  # Input data used by miex
    ├── README.me
    ├── example1.input                           # Exemplary input file
    ├── example2.input                           # Exemplary input file
    ├── miex.py                                  # Source code of miex
    ├── miex_app.py                              # Python script to run miex via streamlit
    ├── miex_example.py                          # Exemplary python script on how to use miex in a user build python code
    ├── miex_notebook.ipynb                      # Jupyter notebook on how to use miex
    ├── requirements.txt                         # Required python packages for miex
    ├── run_miex.py                              # Python script to run miex with input files
    └── test_miex.py                             # Python script for test purposes


## Copyright

The original source code written in `FORTRAN90` is distributed under the [CPC license](https://www.elsevier.com/about/policies/open-access-licenses/elsevier-user-license/cpc-license).
