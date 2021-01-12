import sys
sys.path.append("..")

import numpy as _np
from pyoifits import utils as _u   

import pyoifits
import numpy as np

pyoifits.set_merge_settings(station_distance=1e9, array_distance=1e9)

gravity = pyoifits.open("test1.fits")
pionier = pyoifits.open("test2.fits")
print('Verify GRAVITY & PIONIER')
gravity.verify('fix+warn')
pionier.verify('fix+warn')
merged = gravity + pionier
print('Verify merged file')
merged.verify('fix+warn')



# WGS coordinates of nominal VLTI centre  and nominal station positions.
# (InterfaceControl Document between VLTI and its Instruments (Part I)
# Document ESO-045686, v. 7.3)
arrname = 'VLTI'
lat = -24.62743941
lon = -70.40498689
alt = 2669
sta_name = ['A0', 'B2', 'D0', 'C1']
tel_name = ['AT1', 'AT2', 'AT3', 'AT4']
diameter = [1.8, 1.8, 1.8, 1.8]
staenu = gravity[1].STAXYZ
staenu[:,[0,1]] = -staenu[:,[0,1]]

# Build the OI_ARRAY table 
# Targets
simbad_id = ['beta Car', 'η Car']
category = ['SCI', 'CAL']
array  = pyoifits.new_array_hdu(arrname=arrname, lat=lat, lon=lon, alt=alt,
            tel_name=tel_name, sta_name=sta_name, staenu=staenu,
            diameter=1.8)
target = pyoifits.new_target_hdu_from_simbad(simbad_id, category=category)


# Sample wavelengths
insname = 'TESTING-3CHANNELS'
wave = [2.0e-6, 2.2e-6, 2.4e-6]
band = [0.2e-6, 0.2e-6, 0.2e-6]
wavelength = pyoifits.new_wavelength_hdu(insname=insname, eff_wave=wave, eff_band=band)


# Visibility amplitudes (assume unresolved)
sta_index = [[1, 2], [1, 3], [1,4], [2,3], [2,4], [3, 4]] * 2
target_id = [1] * 6 + [2] * 6
vis2data = [[1, 1, 1]] * 12
visamp = [[1, 1, 1]] * 12
visphi = [[0,0,0]] * 12
vis2err = [[0.05, 0.05, 0.05]] * 12
mjd = [58880] * 6 + [58880.2] * 6



vis2 = pyoifits.new_vis2_hdu(insname=insname, arrname=arrname, mjd=mjd, 
        target_id=target_id, sta_index=sta_index, vis2data=vis2data, 
        vis2err=vis2err)
vis = pyoifits.new_vis_hdu(insname=insname, arrname=arrname, mjd=mjd,
       target_id=target_id, sta_index=sta_index, visamp=visamp,
        visphi=visphi, amptyp='absolute', phityp='absolute') 

sta_index = [[1,2,3], [1,2,4], [1,3,4], [2,3,4]] * 2
target_id = [1] * 4 + [2] * 4
t3phi = [[0,0,0]] * 8
mjd = [58880] * 4 + [58880.2] * 4


t3 = pyoifits.new_t3_hdu(insname=insname, arrname=arrname, mjd=mjd,
    u1coord=0., v1coord=0., u2coord=0., v2coord=0.,
    target_id=target_id, sta_index=sta_index, t3phi=t3phi)

mjd = [58880] * 4 + [58880.2] * 4
sta_index = [1, 2, 3, 4] * 2
fluxdata = [10, 10, 10, 10, 5, 5, 5, 5]
flux = pyoifits.new_flux_hdu(insname=insname, mjd=mjd, fluxdata=fluxdata, 
    target_id=target_id, calibrated=False, arrname=arrname, sta_index=sta_index)

obs = pyoifits.OIFITS2([pyoifits.PrimaryHDU2(), target, array, wavelength, vis, vis2, t3, flux])
 
