import sys
sys.path.append("..")

import pyoifits
import numpy as np

pyoifits.set_merge_settings(station_distance=1e9, array_distance=1e9)

oifits1 = pyoifits.open("test1.fits")
oifits2 = pyoifits.open("test2.fits")
oifits = oifits1 + oifits2
arr = oifits1[1]

tel = list(arr.TEL_NAME)
sta = list(arr.STA_NAME)
index = list(arr.STA_INDEX)
diameter = list(arr.DIAMETER)
xyz = [list(a) for a in arr.STAXYZ]

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
sta_enu = [[-14.642, -55.812, 4.54],
           [ -1.863, -68.334, 4.54],
           [106.648, -39.444, 4.54],
           [ 59.810,  96.706, 4.54]]

# Build the OI_ARRAY table 
# Targets
simbad_id = ['beta Car', 'η Car']
category = ['SCI', 'CAL']
array  = pyoifits.new_array_hdu(arrname=arrname, lat=lat, lon=lon, alt=alt,
            tel_name=tel_name, sta_name=sta_name, sta_enu=sta_enu,
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
        target_id=target_id, sta_index=sta_index, vis2data=vis2data, 
        vis2err=vis2err)

obs = pyoifits.OIFITS2([pyoifits.PrimaryHDU2(), target, array, wavelength, vis2])
