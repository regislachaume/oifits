import sys

sys.path.append("..")
import pyoifits as oifits
from pyoifits.hdu.data import _DataHDU
from pyoifits.hdu.wavelength import _WavelengthHDU

import numpy as _np
from numpy import ma as _ma
from astropy.io import fits as _fits

files = [f"gravity-{i}.fits" for i in range(1, 5)]
hdulists = [oifits.open(f) for f in files]
for hdulist in hdulists: 
    for fhdu in hdulist.get_fluxHDUs(): 
        pass
        if 'FLUX' in fhdu.columns.names: 
            fhdu.rename_columns(FLUX='FLUXDATA')
    hdulist.verify('silentfix+warn')

hdulist = oifits.merge(*hdulists)
hdulist.verify('fix+warn')

hdulist = hdulist.bin_spectral_channels(50)


#fig = hdulist.visualize('spatial_frequency', 'VIS2DATA', target='CO_Ori_A',
#        insname='GRAVITY_SC')
#fig.show()
