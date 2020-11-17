import sys

sys.path.append("..")
import pyoifits as oifits

files = [f"gravity-{i}.fits" for i in range(1, 5)]
hdulists = [oifits.open(f) for f in files]
for hdulist in hdulists: 
    for fhdu in hdulist.get_fluxHDUs(): 
        if 'FLUX' in fhdu.columns.names: 
            fhdu.rename_columns(FLUX='FLUXDATA')
    hdulist.verify('fix+warn')

hdulist = oifits.merge(*hdulists)
hdulist.verify('fix+warn')

fig = hdulist.visualize('spatial_frequency', 'VIS2DATA', target='CO_Ori_A',
        insname='GRAVITY_SC')
fig.show()
