## Installation

Site-wide installation can be performed with ```sudo -H pip3 install oifits```

## Short example
 
```
>>> import oifits
>>> 
>>> hdulist = oifits.open('oifits.fits')
>>> hdulist.verify('check+warn')
WARNING: VerifyWarning: Verification reported errors: [astropy.io.fits.verify]
WARNING: VerifyWarning: HDU 1: [astropy.io.fits.verify]
WARNING: VerifyWarning:     Column 'TEL_NAME' type should be 8-bytes string but is 3-bytes string.  Ignored. [astropy.io.fits.verify]
WARNING: VerifyWarning:     Column 'STA_NAME' type should be 8-bytes string but is 2-bytes string.  Ignored. [astropy.io.fits.verify]
[...]
>>> hdulist
[<PrimaryHDU>, <ArrayHDU (4R)>, <TargetHDU (2R)>, <WavelengthHDU (210SC=210R)>, <WavelengthHDU (5SC=5R)>, <VisHDU (5SC×6R)>, <Vis2HDU (5SC×6R)>, <T3HDU (5SC×4R)>, <FluxHDU (5SC×4R)>, <VisHDU (210SC×6R)>, <Vis2HDU (210SC×6R)>, <T3HDU (210SC×4R)>, <FluxHDU (210SC×4R)>]
>>> tab = hdulist.to_table()
>>> <Table length=9030>
 TARGET  CHANNEL REF_CHANNEL_BITFIELD   EFF_WAVE    EFF_BAND  U1COORD V1COORD ...   error    INSNAME   ARRNAME STA_CONFIG     MJD     INT_TIME
  str8    int64         int64           float32     float32   float64 float64 ...  float64    str10      str4     str8      float64   float64 
-------- ------- -------------------- ----------- ----------- ------- ------- ... --------- ---------- ------- ---------- ----------- --------
CO_Ori_A       1                   -- 2.02369e-06 8.50000e-08   5.488  15.771 ... 0.0041681 GRAVITY_FT    VLTI      C1-D0 57643.37825  273.000
[...]

```
# oifits
