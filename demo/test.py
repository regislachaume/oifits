import sys
sys.path.append("..")

import pyoifits
import numpy as np

pyoifits.set_merge_settings(station_distance=1e9, array_distance=1e9)

oifits1 = pyoifits.open("test1.fits")
oifits2 = pyoifits.open("test2.fits")


oifits = oifits1 + oifits2

# WGS coordinates of nominal VLTI centre  and nominal station positions.
# (InterfaceControl Document between VLTI and its Instruments (Part I)
# Document ESO-045686, v. 7.3)
arrname = 'VLTI'
lat = -24.62743941
lon = -70.40498689
alt = 2669
sta_name = ['A0', 'B1', 'J1', 'J6']
tel_name = ['AT1', 'AT2', 'AT3', 'AT4']
diameter = [1.8, 1.8, 1.8, 1.8]
staenu = oifits1[1].STAXYZ


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
vis2err = [[0.05, 0.05, 0.05]] * 12
mjd = np.linspace(58880, 58880.2, 12)

vis2 = pyoifits.new_vis2_hdu(insname=insname, arrname=arrname, mjd=mjd, 
        ucoord=0., vcoord=0.,
        target_id=target_id, sta_index=sta_index, vis2data=vis2data, 
        vis2err=vis2err)

obs = pyoifits.OIFITS2([pyoifits.PrimaryHDU2(), target, array, wavelength, vis2])
vis2 = oifits1[-4]
# 
for c in 'XYZ':
    oifits1[1].header['ARRAY' + c] = array.header['ARRAY' + c]
    oifits1[1].STAXYZ = array.STAXYZ

#print('no ref')
#uvw = vis2.get_stauvw()
#print('ref')
#uvw_ref = vis2.get_stauvw(refraction=True)

u, v, w = vis2.get_uvw(refraction=True)
u0, v0 = vis2.UCOORD, vis2.VCOORD
print(u-u0, v-v0)

u, v, w = vis2.get_uvw(refraction=False)
u0, v0 = vis2.UCOORD, vis2.VCOORD
print(u-u0, v-v0)


