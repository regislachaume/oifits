# I generated a lot of OIFITS files with aspro spanning most 
# supported instruments at VLTI, CHARA, NPOI, & SUSI. We read them
# to check our compliance-checker.

import sys
sys.path.append("..")

import os
import pyoifits as oifits

templates_dir = 'templates'
templates = sorted(os.listdir(templates_dir))
templates = [os.path.join('templates', t) for t in templates]

for filename in templates:
    print(f"Opening {filename}")
    with oifits.open(filename) as hdulist:
        print(f"Verifying OIFITS")
        hdulist.verify('fix+warn')  # a few may have minor issues
        print(f"Show contents\n{hdulist}")
    print("\n---\n")

print(f"Opening directory")
hdulist = oifits.openlist(templates)
print(f"Verifying OIFITS")
hdulist.verify('fix+warn') # here nothing because openlist will 
                           # silently fix
print(f"Show contents\n{hdulist}")
