import sys
sys.path.append("..")

import pyoifits

file = '/data/lachaume/Cloud/Dropbox/Work/observ/87A/087.C-0227/PIONIER/oifits/2011-05-21-lambda_Sco_calib.fits'

h1 = pyoifits.open(file)
h2 = h1.to_version(2)

