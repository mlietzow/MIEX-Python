# MIEX for Python

[![arXiv](https://img.shields.io/badge/arXiv-astro--ph%2F0406118-red)](https://arxiv.org/abs/astro-ph/0406118)
[![bibcode](https://img.shields.io/badge/bibcode-2004CoPhC.162..113W-blue)](https://ui.adsabs.harvard.edu/abs/2004CoPhC.162..113W)
[![doi](https://img.shields.io/badge/doi-10.1016%2Fj.cpc.2004.06.070-yellow)](https://doi.org/10.1016/j.cpc.2004.06.070)

is a Mie scattering code for large grains.


### Run miex

[run_miex.py](run_miex.py) executes miex with the parameters of the input file:

```bash
python3 run_miex.py example1.input
python3 run_miex.py example2.input
```
The results are stored in the `results` directory.


### Run miex via streamlit

[miex_app.py](miex_app.py) can be used as a simple Streamlit web app. Install streamlit
```bash
pip3 install streamlit
```
and start the web app:
```bash
streamlit run miex_app.py
```


### Import miex

[miex.py](miex.py) can be imported and used in any python script (see [example script](miex_example.py) or [jupyter notebook](miex_notebook.ipynb)):

```bash
python3 miex_example.py
```

To calculate the efficiency factors and scattering amplitude functions (optionally), use e.g.,
```python
miex.shexqnn2(x=1.0, m=complex(1.5, 0.0), nang=91, doSA=True)
```
where `x` is the size parameter, `m` is the complex refractive index, `nang` is the half number of scattering angles in the intervall $0$ to $\pi/2$, and `doSA` is either `True` (calculate scattering amplitude functions) or `False`.
The function returns the extinction efficiency $Q_{\rm ext}$, absorption efficiency $Q_{\rm abs}$, scattering efficiency $Q_{\rm sca}$, backscattering effiency $Q_{\rm bk}$, radiation pressure effiency $Q_{\rm pr}$, single scattering albedo $A$, scattering assymetry factor $g$, and scattering amplitude functions $S_1$ and $S_2$ (in this order).
The scattering amplitude functions are an array with size `2*nang-1`.

To calculate the scattering matrix elements, use
```python
miex.scattering_matrix_elements(S1, S2)
```
where $S_1$ and $S_2$ are the scattering amplitude functions.
The function returns the scattering matrix elements $S_{11}$, $S_{12}$, $S_{33}$, $S_{34}$ (in this order).


### Test miex

[test_miex.py](test_miex.py) includes some test routines. The results are compared with results by [Bohren & Huffman (1998)](https://doi.org/10.1002/9783527618156) and by [Wiscombe (1979)](https://doi.org/10.5065/D6ZP4414):

```bash
python3 test_miex.py
```
