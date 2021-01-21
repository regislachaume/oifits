"""
A package to read, write, manipulate, and create files containing the Optical 
Interferometry Data Exchange standard format defined by
* OIFITS 2: Gilles Duvert, John Young & Christian Hummel (2017),
  Astronomy & Astrophysics, 597, A8, doi:10.1051/0004-6361/201526405, 
  arXiv:1510.04556 
* OIFITS 1: T.A. Pauls et al. (2005), Publications of the Astronomical 
  Society of the Pacific 117, 1255, doi:10.1086/444523, 
  arXiv:astro-ph/050818

Typical use

>>> from pyoifits import open as oifits_open
>>> hdulist = oifits_open('filename.fits')

"""

from .oifits import *

__version__ = "0.4.7"
__author__ = "Régis Lachaume"
__email__ = "regis.lachaume@gmail.com"
__date__ = "2020-01-21"
