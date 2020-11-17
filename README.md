## Purpose

Read and manipulate Optical Interferometry FITS files.  For a definition of
the standards 
* version 2: Duvert at al. (2017), A&A 597, A8 ([abstract](https://ui.adsabs.harvard.edu/abs/2017A%26A...597A...8D/abstract "ADS link")) ([pre-print](https://arxiv.org/pdf/1510.04556 "PDF on arxiv"))
* version 1: Pauls et al. (2005), PASP 117, 1255 ([abstract](https://ui.adsabs.harvard.edu/abs/2005PASP..117.1255P/abstract "ADS link")) ([pre-print](https://arxiv.org/pdf/astro-ph/0508185 "PDF on arxiv"))

## Installation

Site-wide installation will be performed with `sudo -H pip3 install pyoifits` on unix-like systems.

At a user level, within a [virtual environment](https://docs.python.org/3/library/venv.html "venv package"), `pip3 install pyoifits`. 

## Short example

Read and merge to OIFITS datasets and tranform to a standard [`astropy`](https://www.astropy.org) table.
    
```python

import oifits

data1 = oifits.read('file1.fits')
data2 = oifits.read('file2.fits')

data = data1 + data2

tab, C = data.to_table(correlations='numpy')
val = tab['value']
err = tab['error']
covar = C * np.outer(err, err)
```

There is also a short [demo](https://github.com/loqueelvientoajuarez/oifits/blob/master/demo/intro.ipynb "Jupyter notebook demo").

## Release notes

### 0.4

New features:
* Method `to_table()` optionally returns a correlation matrix

### 0.4.1

Bug fixes:
* `WAVELMIN` & `WAVELMAX` FITS keywords have a correct number of significant digits
* `EQUINOX` FITS keyword now passes verification
* No runtime error when a duplicate `INSNAME` has to be renamed
* OI data (`OI_VIS`, `OI_VIS2`, `OI_T3`, `OI_FLUX`) tables merge correctly
* `OI_VIS` tables with different `PHITYP`, `AMPTYP`, `AMPORDER`, or `PHIORDER` can no longer be merged
* `OI_FLUX` tables with different `FOV` or  `FOVTYPE` can not longer be merged

New features:
* `OIFITS1` & `OIFITS2` classes have `get_visHDUs`, `get_vis2HDUs`, `get_T3HDUs` methods

### 0.4.2

New features:
* `to_table()` allows the user to restrict to a subset of target names, 
instrument setup names, arrays, wavelength range, and/or date range
* OIFITS classes have new method `visualize()` to plot an interferometric
observable as a function of time, baseline length, or spatial frequency.
