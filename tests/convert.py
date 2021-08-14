#! /usr/bin/env python3

import sys
import os

sys.path.append("..")

import pyoifits
import numpy as np

v1 = os.path.join('fits', 'convert-v1.fits')
v2 = os.path.join('fits', 'convert-v2.fits')

h1 = pyoifits.open(v1)
# mess a bit with invalid errors and TARGET_ID
h1[5].T3AMPERR[:,-1] = np.nan
h1[5].T3AMPERR[:,0] = -1
h1[1]._substitute_ids({8: -1})

h2 = h1.to_version(2)
h2.writeto(v2, overwrite=True)
