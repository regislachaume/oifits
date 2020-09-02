import sys
sys.path.append("..")

import pyoifits

pyoifits.set_merge_settings(station_distance=1e9, array_distance=1e9)

oifits1 = pyoifits.open("test1.fits")
oifits2 = pyoifits.open("test2.fits")
oifits = oifits1 + oifits2
arr = oifits1[1]
print(type(arr.header))

tel = list(arr.TEL_NAME)
sta = list(arr.STA_NAME)
index = list(arr.STA_INDEX)
diameter = list(arr.DIAMETER)
xyz = [list(a) for a in arr.STAXYZ]

frame = arr.header['FRAME']
x, y, z = [arr.header[n] for n in ['ARRAYX', 'ARRAYY', 'ARRAYZ']]
arrname = arr.header['ARRNAME']

print('build new table')
tab  = pyoifits.new_array_hdu('VLTI',
            lat=-24.62743941, lon=-70.40498689, alt=2669, frame=frame,
            tel_name=tel, sta_name=sta, sta_enu=xyz)

