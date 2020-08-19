## Purpose

Read and manipulate Optical Interferometry FITS files.  For a definition of
the standards 
* version 2: Duvert at al. (2017), A&A 597, A8 [abstract](https://ui.adsabs.harvard.edu/abs/2017A%26A...597A...8D/abstract "ADS link") [pre-print](https://arxiv.org/pdf/1510.04556 "PDF on arxiv")
* version 1: Pauls et al. (2005), PASP 117, 1255 [abstract](https://ui.adsabs.harvard.edu/abs/2005PASP..117.1255P/abstract "ADS link") [pre-print](https://arxiv.org/pdf/astro-ph/0508185 "PDF on arxiv")

## Installation

Site-wide installation will be performed with `sudo -H pip3 install oifits`.

At a user level, within a virtual environment, `pip3 install oifits`. 

## Short example

Read and merge to OIFITS datasets and tranform to a standard ``astropy`` table.
    
```python

import oifits

data1 = oifits.read('file1.fits')
data2 = oifits.read('file2.fits')

data = data1 + data2

tab = data.to_table()
```

[Demo on github](https://github.com/loqueelvientoajuarez/oifits/blob/master/demo/intro.ipynb Jupyter notebook demo)
