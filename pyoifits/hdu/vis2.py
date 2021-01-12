from .data import _DataHDU11, _DataHDU22
from .t2 import _T2HDU
from .. import utils as _u
from .wavelength import _NW

import numpy as _np

class _Vis2HDU(_T2HDU):
    """

OI_VIS2 binary table extension containing the squared visibility amplitudes

    """

    _EXTNAME = 'OI_VIS2'
    _COLUMNS = [
        ('VIS2DATA', True, 'D', (_NW,), None, None, None,
            'squared visibility amplitude'), 
        ('VIS2ERR',  True, 'D', (_NW,), None, None, None,
            'uncertaintly on squared visibility amplitude'),
    ]
    
    def get_obs_type(self, name, shape='data', flatten=False):
        
        return self._resize_data('absolute', shape, flatten)

    @classmethod
    def from_data(cls, *, insname, mjd, vis2data, target_id, sta_index,
            fits_keywords={}, **columns):
        """

Build an OI_VIS2 table from data.  In the following NWAVE indicates the
number of spectral channels in the configuration and NOBS the number
of observations.

Arguments
---------

version (int)
    version of the OIFITS standard if it cannot be determined from
    context (optional)

date (str)
    date at the start of the first observation (optional, deduced from
    mdj and int_time)
insname (str)
    name of the OI_WAVELENGTH table containing the instrumental setup
arrname (str)         
    name of the OI_ARRAY table containing array properties (optional in
    version 1)
fits_keywords (dict)
    additional FITS keywords (optional)

mjd (float, NOBS)
    MJD of observation
int_time (float, NOBS)
    Integration time
target_id (int, NOBS)
    Target ID in the OI_TARGET table
sta_index (int, NOBS × 2)
    Station IDs of the telescope pair in the OI_ARRAY target (optional in
    version 1)
ucoord, vcoord (float, NOBS)
    (u, v) coordinates of the  baselines

vis2data (float, NOBS × NWAVE)
    squared visibility amplitudes
vis2err (float, NOBS × NWAVE)
    uncertainty on the latter (optional, defaults to 0.)
flag (bool, NOBS × NWAVE)
    flag for data (optional, defaults to False)

Arguments in OIFITS1 only
-------------------------

time (float, NOBS)
    Seconds since midnight UT on first observation (DATE-OBS i.e. date)

Arguments in OIFITS2 only
-------------------------

corrindx_vis2data (int, NOBS × NWAVE)
    index in the correlation matrix for the first vis2data (optional)

Additional arguments
--------------------

Any additional keyword argument will be appended as a non-standard FITS 
column with its name prefixed with NS_ 

        """
        columns = dict(vis2data=vis2data, target_id=target_id,
                sta_index=sta_index, **columns)

        return super().from_data(insname=insname, mjd=mjd,
                    fits_keywords=fits_keywords, **columns)

class Vis2HDU1(
        _Vis2HDU,
        _DataHDU11,
      ):
    """

First revision of the OI_VIS2 binary table, OIFITS v. 1

    """
    pass

class Vis2HDU2(
        _Vis2HDU,
        _DataHDU22,
      ):
    """ 

Second revision of the OI_VIS2 binary table, OIFITS v. 2

    """
    _COLUMNS = [
        ('CORRINDX_VIS2DATA', False, '1J', (), _u.is_strictpos, None, None,
            'index of 1st visib. in matching OI_CORR matrix'),
    ]

new_vis2_hdu = _Vis2HDU.from_data
