## Purpose

Read and manipulate Optical Interferometry FITS files.  See Duvert at al.
(2007, A&A 597, A8) for a full definition of the standard. 

## Installation

Site-wide installation will be performed with ```sudo -H pip3 install oifits``` (when this package is ready)

## Short example

Read and merge to OIFITS datasets and tranform to a standard ```astropy``` table.

```import oifits

data1 = oifits.read('file1.fits')
data2 = oifits.read('file2.fits')
data = data1 + data2
tab = data.to_table()
```
