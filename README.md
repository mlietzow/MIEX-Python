# MIEX-Python

is a Mie scattering code for large grains written in Python and based on [MIEX](https://ui.adsabs.harvard.edu/abs/2018ascl.soft10019W) by [Wolf & Voshchinnikov (2004)](https://ui.adsabs.harvard.edu/abs/2004CoPhC.162..113W).

The following quantities for

1. single grain sizes / chemical components and
2. mixtures of chemically different grains with a size distribution

can be calculated:

- Scattering matrix elements $S_{11}$, $S_{12}$, $S_{33}$, and $S_{34}$,
- Extinction efficiency factor ($Q_\mathrm{ext}$) and Extinction cross-section ($C_\mathrm{ext}$),
- Scattering efficiency factor ($Q_\mathrm{sca}$) and Scattering cross-section ($C_\mathrm{sca}$),
- Absorption efficiency factor ($Q_\mathrm{abs}$) and Absorption cross-section ($C_\mathrm{abs}$),
- Backscattering efficiency factor ($Q_\mathrm{bk}$) and Backscattering cross-section ($C_\mathrm{bk}$),
- Radiation pressure efficiency factor ($Q_\mathrm{pr}$),
- Albedo,
- Scattering asymmetry factor ($g$).

Modified and ported to Python with permission of S. Wolf.
For the original source code of MIEX written in `FORTRAN90`, we refer to

[![arXiv](https://img.shields.io/badge/arXiv-astro--ph%2F0406118-b31b1b)](https://arxiv.org/abs/astro-ph/0406118)
[![ascl](https://img.shields.io/badge/ascl-1810.019-262255)](https://ascl.net/1810.019)
[![bibcode](https://img.shields.io/badge/bibcode-2004CoPhC.162..113W-1c459b)](https://ui.adsabs.harvard.edu/abs/2004CoPhC.162..113W)
[![doi](https://img.shields.io/badge/doi-10.1016%2Fj.cpc.2004.06.070-fab70c)](https://doi.org/10.1016/j.cpc.2004.06.070)


## Requirements

To run MIEX-Python, at least **Python 3.6** and the following packages are required:
 - numba
 - numpy

To check the current version or to install the packages, run:

```bash
python3 --version
pip3 install -r requirements.txt
```

## Installation

The package can be installed via `pip`, change into the directory `MIEX-Python` and run:
```
pip install .
```


## Use MIEX-Python

### Run MIEX-Python via streamlit

[miex_app.py](miex_app.py) can be used as a simple Streamlit web app: https://miex-python.streamlit.app/.

Alternatively, install Streamlit

```bash
pip3 install streamlit
```

and start the web app locally:

```bash
streamlit run miex_app.py
```


### Import MIEX-Python
After installation, MIEX-Python can be imported via `import miex` and used in any python script (see [jupyter notebook 1](miex_notebook_1.ipynb)).

To calculate the efficiency factors and scattering amplitude functions (optionally), use e.g.,

```python
result = miex.get_mie_coefficients(x=1.0, m=complex(1.5, 0.0), nang=91, doSA=True)
```

| Variable      | Input Parameter                                          | Type            |
| ------------- | -------------------------------------------------------- | --------------- |
| `x`           | Size parameter                                           | float           |
| `m`           | Complex refractive index                                 | complex         |
| `nang=2`      | Half number of scattering angles in the interval [0, 90] | int, optional   |
| `doSA=False`  | Calculate scattering amplitude functions                 | bool, optional  |
| `nterm=2e7`   | Maximum number of terms to be considered                 | int, optional   |
| `eps=1.0e-20` | Accuracy to be achieved                                  | float, optional |
| `xmin=1.0e-6` | Minimum size parameter                                   | float, optional |

The function returns a dictionary with the following entries:

| Entry         | Output Parameter              | Type                |
| ------------- | ------------------------------| ------------------- |
| `Q_ext`       | Extinction efficiency         | float               |
| `Q_abs`       | Absorption efficiency         | float               |
| `Q_sca`       | Scattering efficiency         | float               |
| `Q_bk`        | Backscattering efficiency     | float               |
| `Q_pr`        | Radiation pressure efficiency | float               |
| `Albedo`      | Single scattering albedo      | float               |
| `g_sca`       | Scattering asymmetry factor   | float               |
| `SA_1`        | Scattering amplitude function | ndarray, complex    |
| `SA_2`        | Scattering amplitude function | ndarray, complex    |
| `theta`       | Scattering angles [rad]       | ndarray             |

The scattering amplitude functions are an array with size `2*nang-1`.

To calculate the scattering matrix elements, use

```python
scat_mat = miex.get_scattering_matrix_elements(S1, S2)
```

where $S_1$ and $S_2$​ are the scattering amplitude functions.
The function returns a dictionary with the scattering matrix elements $S_{11}$​, $S_{12}$​, $S_{33}$​, $S_{34}$ as entries for instance `scat_mat["S_11"]`.


### Test MIEX-Python

[test_miex.py](test_miex.py) includes some test routines. The results are compared with results by [Bohren & Huffman (1998)](https://doi.org/10.1002/9783527618156) and by [Wiscombe (1979)](https://doi.org/10.5065/D6ZP4414):

```bash
python3 test_miex.py
```


## Project structure

    .
    ├── ri-data                                  # Input data used by MIEX-Python
    ├── README.me
    ├── miex
    │   └── miex.py                              # Source code of MIEX-Python
    ├── miex_app.py                              # Python script to run MIEX-Python via Streamlit
    ├── miex_notebook_1.ipynb                    # Jupyter notebook on how to use MIEX-Python
    ├── requirements.txt                         # Required python packages for MIEX-Python
    └── test_miex.py                             # Python script for test purposes


## Copyright

The original source code of MIEX written in `FORTRAN90` is distributed under the [CPC license](https://www.elsevier.com/about/policies/open-access-licenses/elsevier-user-license/cpc-license).
Modified and ported to Python with permission of S. Wolf.
