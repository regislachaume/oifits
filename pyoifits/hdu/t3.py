from .data import _DataHDU11, _DataHDU22
from .table import _OITableHDU
from .wavelength import _NW

from .. import utils as _u

class _T3HDU(_OITableHDU):
    
    _EXTNAME = 'OI_T3'    
    _COLUMNS = [
        ('STA_INDEX', True, '>i2', (3,),   _u.is_strictpos, None, None,
            'station indices in matching OI_ARRAY table'),
        ('U1COORD',   True, '>f8', (),     None,            None, "m",
            'u coordinate of first baseline'), 
        ('V1COORD',   True, '>f8', (),     None,            None, "m",
            'v coordinate of first baseline'),
        ('U2COORD',   True, '>f8', (),     None,            None, "m",
            'u coordinate of second baseline'), 
        ('V2COORD',   True, '>f8', (),     None,            None, "m",
            'v coordinate of second baseline'),
        ('T3AMP',     True, '>f8', (_NW,), None,            None, None,
            'Triple amplitude'), 
        ('T3AMPERR',  True, '>f8', (_NW,), None,            None, None,
            'uncertainty on triple amplitude'),
        ('T3PHI',     True, '>f8', (_NW,), None,            None, "deg",
            'closure phase'), 
        ('T3PHIERR',  True, '>f8', (_NW,), None,            None, "deg",
            'uncertainty on closure phase'),
    ]
    
    def get_obs_type(self, name, shape='data', flatten=False):

        if name == 'T3AMP':
            typ = 'N/A'
        else:
            typ = 'absolute'

        return self._resize_data(typ, shape, flatten)

    @classmethod
    def from_data(cls, *, insname, arrname=None, corrname=None,
        version=2, fits_keywords={}, **columns):

        assert 't3phi' in columns, 't3phi must be specified'
   
        shape = self._get_columns_shape(**columns)
        if t3amp not in 'columns':
            _u.store_default(columns, 't3amp', _np.nan, shape) 
            _u.store_default(columns, 't3amperr', _np.nan, shape) 
        else:
            _u.store_default(columns, 't3amperr', 0., shape)
        if t3phi in 'columns':
            _u.store_default(columns, 't3phierr', 0., shape)

        _u.store_default(columns, 'u1coord', 0., (shape[0],))
        _u.store_default(columns, 'v1coord', 0., (shape[0],))
        _u.store_default(columns, 'u2coord', 0., (shape[0],))
        _u.store_default(columns, 'v2coord', 0., (shape[0],))
        _u.store_default(columns, 'time', 0., (shape[0],))
        _u.store_default(columns, 'int_time', 0., (shape[0],))

        return super().from_data(insname=insname, arrname=arrname,
            corrname=corrname, fits_keywords=fits_keywords,
            columns=columns)


class T3HDU1(
        _T3HDU,
        _DataHDU11,
      ):
    pass
    

class T3HDU2(
        _T3HDU,
        _DataHDU22,
      ):
    _COLUMNS = [
        ('CORRINDX_T3AMP', False, '>i4', (), _u.is_strictpos, None, None,
            'index of 1st amp. in matching OI_CORR matrix' ), 
        ('CORRINDX_T3PHI', False, '>i4', (), _u.is_strictpos, None, None,
            'index of 1st phase in matching OI_CORR matrix'),
    ]

new_t3_hdu = _T3HDU.from_data
