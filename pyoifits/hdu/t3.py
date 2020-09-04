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
