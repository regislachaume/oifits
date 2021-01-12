from .data import _DataHDU, _DataHDU11, _DataHDU22
from .wavelength import _NW

import numpy as _np

from .. import utils as _u

class _T3HDU(_DataHDU):
    """
OI_T3 binary table extension containing closure phases and triple 
product amplitudes.
    """
    
    _EXTNAME = 'OI_T3'    
    _COLUMNS = [
        ('STA_INDEX', True, '3I', (3,),   _u.is_strictpos, None, None,
            'station indices in matching OI_ARRAY table'),
        ('U1COORD',   True, '1D', (),     None,            None, "m",
            'u coordinate of first baseline'), 
        ('V1COORD',   True, '1D', (),     None,            None, "m",
            'v coordinate of first baseline'),
        ('U2COORD',   True, '1D', (),     None,            None, "m",
            'u coordinate of second baseline'), 
        ('V2COORD',   True, '1D', (),     None,            None, "m",
            'v coordinate of second baseline'),
        ('T3AMP',     True, 'D',  (_NW,), None,            None, None,
            'Triple amplitude'), 
        ('T3AMPERR',  True, 'D',  (_NW,), None,            None, None,
            'uncertainty on triple amplitude'),
        ('T3PHI',     True, 'D',  (_NW,), None,            None, "deg",
            'closure phase'), 
        ('T3PHIERR',  True, 'D',  (_NW,), None,            None, "deg",
            'uncertainty on closure phase'),
    ]
    
    def get_obs_type(self, name, shape='data', flatten=False):

        type = 'absolute'

        return self._resize_data(type, shape, flatten)

    def get_uvw(self, refraction=False):
        """

Get the (u, v, w) coordinates of the first two baselines

Arguments
---------

refraction (bool, optional, default: False)
    Whether refraction is taken into account

Returns
-------

(u1, v1, w1), (u2, v2, w2)
    Projected baselines in metres

Precision
---------
    A few centimetres for hectometric baselines, mainly due to the uncertainty
    on the time-averaging.
    

        """
        uvw = self.get_stauvw(refraction=refraction)
        uvw1 = uvw[:,1,:] - uvw[:,0,:]
        uvw2 = uvw[:,2,:] - uvw[:,1,:]
      
        return uvw1.T, uvw2.T

    def update_uv(self):

        if self.get_arrayHDU() is None:
            return

        (u1, v1, w1), (u2, v2, w2) = self.get_uvw()
        self.U1COORD = u1
        self.V1COORD = v1
        self.U2COORD = u2
        self.V2COORD = v2
        
    @classmethod
    def from_data(cls, *, insname, mjd, target_id, sta_index, 
            t3phi, fits_keywords={}, **columns):
        """

Build an OI_T3 table from data.  In the following NWAVE indicates the
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
sta_index (int, NOBS × 3)
    Station IDs of the telescope triplet in the OI_ARRAY target (optional 
    in version 1)
u1coord, v1coord, u2coord, v2coord (float, NOBS)
    (u, v) coordinates of the first two baselines.
t3phi (float, NOBS × NWAVE)
    closure phase
t3phierr (float, NOBS × NWAVE)
    uncertainty on the latter (optional, defaults to 0.)
t3amp (float, NOBS × NWAVE)
    triple product amplitude (optional, defaults to NULL) 
t3phierr (float, NOBS × NWAVE)
    uncertainty on the latter (optional, defaults to 0.)
flag (bool, NOBS × NWAVE)
    flag for data (optional, defaults to False)

Arguments in OIFITS1 only
-------------------------

time (float, NOBS)
     Seconds since midnight UT on first observation (DATE-OBS i.e. date)


Arguments in OIFITS2 only
-------------------------

corrindx_t3phi (int, NOBS × NWAVE)
    index in the correlation matrix for the first t3phi (optional)
corrindx_t3amp (int, NOBS × NWAVE)
    index in the correlation matrix for the first t3amp (optional)

Additional arguments
--------------------

Any additional keyword argument will be appended as a non-standard FITS 
column with its name prefixed with NS_ 

        """
        t3amp = _np.empty_like(t3phi)
        t3amp[...] = _np.nan
        _u.store_default(columns, 't3amp', default=t3amp)
        _u.store_default(columns, 'u1coord', default=0.)
        _u.store_default(columns, 'u2coord', default=0.)
        _u.store_default(columns, 'v1coord', default=0.)
        _u.store_default(columns, 'v2coord', default=0.)

        columns = dict(t3phi=t3phi, target_id=target_id, 
                    sta_index=sta_index, **columns)

        return super().from_data(insname=insname, mjd=mjd,
                fits_keywords=fits_keywords, **columns)

class T3HDU1(
        _T3HDU,
        _DataHDU11,
      ):
    """

First revision of the OI_T3 binary table, OIFITS v. 1

    """
    pass
    

class T3HDU2(
        _T3HDU,
        _DataHDU22,
      ):
    """

Second revision of the OI_T3 binary table, OIFITS v. 2

    """
    _COLUMNS = [
        ('CORRINDX_T3AMP', False, '1J', (), _u.is_strictpos, None, None,
            'index of 1st amp. in matching OI_CORR matrix' ), 
        ('CORRINDX_T3PHI', False, '1J', (), _u.is_strictpos, None, None,
            'index of 1st phase in matching OI_CORR matrix'),
    ]

    pass


new_t3_hdu = _T3HDU.from_data
