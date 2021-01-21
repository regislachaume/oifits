## Brief intro

### Purpose

Read and manipulate Optical Interferometry FITS files version 1 and 2.  

For a definition of the standards 
* version 2: Duvert at al. (2017), A&A 597, A8 ([abstract](https://ui.adsabs.harvard.edu/abs/2017A%26A...597A...8D/abstract "ADS link")) ([pre-print](https://arxiv.org/pdf/1510.04556 "PDF on arxiv"))
* version 1: Pauls et al. (2005), PASP 117, 1255 ([abstract](https://ui.adsabs.harvard.edu/abs/2005PASP..117.1255P/abstract "ADS link")) ([pre-print](https://arxiv.org/pdf/astro-ph/0508185 "PDF on arxiv"))

### Installation

It is packaged as [pyoifits](https://pypi.org/project/pyoifits/).

Site-wide installation will be performed with `sudo -H pip3 install pyoifits` on unix-like systems. At a user level, within a [virtual environment](https://docs.python.org/3/library/venv.html "venv package"), `pip3 install pyoifits`. 

### Short example

Read and merge to OIFITS datasets and tranform to a standard [`astropy`](https://www.astropy.org) table.
    
```python

import pyoifits as oifits

data = oifits.openlist(['file1.fits', 'file2.fits'])

tab, corr = data.to_table(correlations='numpy')
obs = tab['observable']
x = tab['value']
dx = tab['error']
cov = corr * np.outer(dx, dx)
```

There is also a short [demo](https://github.com/loqueelvientoajuarez/oifits/blob/master/demo/intro.ipynb "Jupyter notebook demo").

## Classes and functions

### Classes

Classes are in `astropy.io.fits` tree, they derive either from `HDUList` (FITS file) or `BinTableHDU` (FITS table).  The end user shouldn't need to manipulate them explicitly as helper functions are provided.


An `OIFITS1` object (OIFITS standard version 1, < `_OIFITS` < `astropy.io.fits.HDUList`) contains
1. one `PrimaryHDU1` object (< `astropy.io.fits.PrimaryHDU`)
2. several `_OITableHDU1` objects (abstract OI table class, < `_OITableHDU` < `astropy.io.fits.BinTableHDU`) which can be of the derived types
    * `TargetHDU1` (`OI_TARGET` extrension) 
    * `ArrayHDU1` (`OI_ARRAY` extension)
    * `WavelengthHDU1` (`OI_WAVELENGTH` extension, < `astropy.io.fits.BinTableHDU`)
    * `_DataHDU1` (abstract data class, < `astropy.io.fits.BinTableHDU`) 
        * `VisHDU1` (`OI_VIS` extension)
        * `Vis2HDU1` (`OI_VIS2` extension)
        * `T3HDU1` (`OI_T3` extension)
3. any other FITS extension 

An `OIFITS2` object (OIFITS standard version 2, < `_OIFITS` < `astropy.io.fits.HDUList`) contains
1. one `PrimaryHDU2` object (< `astropy.io.fits.PrimaryHDU`)
2. several `_OITableHDU2` (abstract OI table class, < `_OITableHDU` < `astropy.io.fits.BinTableHDU`) which can be of the derived types
    * `TargetHDU2` (`OI_TARGET` extension) 
    * `ArrayHDU2` (`OI_ARRAY` extension)
    * `WavelengthHDU2` (`OI_WAVELENGTH` extension)
    * `_DataHDU2` (abstract data class)
        * `VisHDU2` (`OI_VIS` extension)
        * `Vis2HDU2` (`OI_VIS2` extension)
        * `T3HDU2` (`OI_T3` extension)
        * `FluxHDU1` (`OI_FLUX` extension)
    * `CorrHDU1` (`OI_CORR` extension)
    * `InspolHDU1` (`OI_INSPOL` extension)
3. any other FITS extension

### Functions

* Reading and manipulating OIFITS
    * `open` (open an OIFITS file)
    * `openlist` (open a list of files and merge them)
    * `merge` (merge several OIFITS)
    * `set_merge_settings` (determine how duplicate targets/stations are merged)
* Creating OI tables from scratch
    * `new_target_hdu` (create an `OI_TARGET` extension)
    * `new_target_hdu_from_simbad`
    * `new_array_hdu` (create an `OI_ARRAY` extension)
    * `new_wavelength_hdu` (create an `OI_WAVELENGTH` extension)
    * `new_vis_hdu` (create an `OI_VIS` extension)
    * `new_vis2_hdu` (create an `OI_VIS2` extension)
    * `new_t3_hdu` (create an `OI_T3` extension)
    * `new_flux_hdu` (create an `OI_FLUX` extension)

### Methods

* `_OIFITS` (`OIFITS1` and `OIFITS2`)
    * manipulation
        * `bin_spectral_channels` (downgrade spectral resolution)
        * `trim` (keep only wavelengths, targets, instruments, ... of interest)
        * `to_table` (transform to a table with one scalar observable per line)
        * `to_version` (transform between versions of the OIFITS standard)
    * visualisation
        * `visualize` (quick plot)
    * data updating
        * `verify` (check standard compliance and mend if possible)
        * `update_primary_header` (update primary header using the info in other tables)
        * `update_uv` (compute spatial frequencies using array data)
    * contents listing 
        * `get_OITableHDUs`
        * `get_targetHDU`
        * `get_arrayHDUs` & `get_arrayHDU`
        * `get_wavelengthHDUs` & `get_wavelengthHDU`
        * `get_dataHDUs` 
        * `get_visHDUs` 
        * `get_vis2HDUs`
        * `get_t3HDUs`
        * (`OIFITS2` only) `get_corrHDUs` & `get_corrHDU`
        * (`OIFITS2` only) `get_fluxHDUs`
        * (`OIFITS2` only) `get_inspolHDUs` & `get_inspolHDU`
* `_DataHDU` (all classes containing OI data)
    * target information
        * `get_target`
        * `get_dec`
        * `get_ra`
        * `get_equinox`
        * `get_parallax`
        * `get_pmra`
        * `get_pmdec`
        * `get_rv`
        * `get_skycoord`
        * `get_spectype`
    * wavelength and spatial frequency information
        * `get_wave`
        * `get_band`
        * `get_uv`
        * `get_uvw`
    * array information
        * `get_location`
        * `get_sta_name`
        * `get_staenu`
        * `get_stauvw`
        * `get_staxyz`  
        * `get_sta_config`
        * `get_tel_name`
        * `get_tel_config`
    * polarisation information
        * (`_DataHDU2` only) `get_jones_matrix`

### Direct data access

Each column of the OI FITS standard can be directly accessed via its standard
name, for instance, `h.UCOORD` or `h.VIS2DATA` are synonyms to `h.data['UCOORD']` and `h.data['VIS2DATA']`, respectively.

## Related projects

* [oifits](https://pypi.org/project/oifits/) reads an OIFITS file into a table, geared towards the Event Horizon Telescope
* [oifits](https://github.com/pboley/oifits) reads and write OIFITS v. 1, mostly geared towards VLTI/MATISSE
* [pyhdust](https://pypi.org/project/pyhdust/) has a read/write module for OIFITS version 1

## Release notes

### 0.4

New features:
* Method `to_table()` optionally returns a correlation matrix
GG
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
* Method `to_table()` allows the user to restrict to a subset of target names, 
instrument setup names, arrays, wavelength range, and/or date range
* OIFITS classes have new method `visualize()` to plot an interferometric
observable as a function of time, baseline length, or spatial frequency.

### 0.4.3

Bug fixes:
* Successive column renamings no longer lead to `KeyError`
* fix `get_uv` 
* Typos in this doc

New features:
* Method `bin_spectral_channels(R)` allows the user to rebin to a given spectral resolution.

### 0.4.4
* Documentation
* Fixed numerous bugs in `get_jones_matrix` and `InspolHDU1`

### 0.4.5

Bug fixes
* Fixed error in creation of `FluxHDU` and `T3HDU` using `from_data`
* Fixed error in creation of OI tables with non-standard columns using `from_data`
* Typos in documentation

New features
* New function `trim()` to keep only wavelengths, instruments, instrumental setups, arrays, targets of interest

### 0.4.6.

Bug fixes
* Fixed `TDIM`=1 keyword error when `trim()` yields N×1 data
* Fixed an issue with `astroquery` returning `str` instead of `bytes`

### 0.4.7

Bug fixes
* Fixed float format for `STAXYZ` column in `ArrayHDU` extension
* `verify()` can now fix strictly negative errors
* duplicate fits warnings in `verify()` are now printed

Testing
* A quick test of `verify()` on a set of JMMC/aspro-generated OIFITS
