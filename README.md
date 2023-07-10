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


### Test miex

[test_miex.py](test_miex.py) includes some test routines. The results are compared with results by [Bohren & Huffman (1998)](https://doi.org/10.1002/9783527618156) and by [Wiscombe (1979)](https://doi.org/10.5065/D6ZP4414):

```bash
python3 test_miex.py
```
