# propitiatory invocation (i.e. the user hasn't installed pyoifits as a package)
import sys
sys.path.append('..') 

import numpy as np
from astropy.table import Table

import pyoifits as oifits

# OIFITS version 2 from GRAVITY @ VLTI
gravity = oifits.open('introfiles/gravity-1.fits')
# pionier = oifits.open('introfiles/pionier-1.fits')
