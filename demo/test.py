import sys
sys.path.append('..')

import oifits

# Open a FITS file produced by the GRAVITY pipeline
hdulist = oifits.oifitsopen('test1.fits')
for h in hdulist[8::4]:
    h.rename_columns(FLUX='FLUXDATA')

# Check compliance
hdulist.verify('fix+warn')

# Transform into a flat table (one scalar observable per line).
tab = hdulist.to_table()
print(tab)
