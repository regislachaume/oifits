#! /usr/bin/env python3

import sys
import os

sys.path.append("..")

import pyoifits
import numpy as np

v1 = os.path.join('fits', 'convert-v1.fits')
v2 = os.path.join('fits', 'convert-v2.fits')

h1 = pyoifits.open(v1)

# mess a bit with invalid errors 
h1[5].T3AMPERR[:,-1] = np.nan
h1[5].T3AMPERR[:,0] = -1

# give invalid version 2 TARGET_ID / STA_ID (<= 0)
thdu = h1.get_targetHDU()
for h in [thdu, *thdu.get_referrers()]:
    h.TARGET_ID[h.TARGET_ID == 8] = -1

sub = {1: 0, 2: -1}
ahdu = h1.get_arrayHDUs()[0]
for h in [ahdu, *ahdu.get_data_referrers()]:
    for old, new in sub.items():   
        h.STA_INDEX[h.STA_INDEX == old] = new

#h2 = h1.to_version(2)
#h2.writeto(v2, overwrite=True)
