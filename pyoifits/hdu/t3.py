from .data import _DataHDU11, _DataHDU22
from .table import _OITableHDU
from .wavelength import _NW

from .. import utils as _u

__all__ = ["T3HDU1", "T3HDU2"]

class _T3HDU(_OITableHDU):
    
    _EXTNAME = 'OI_T3'    
    _COLUMNS = [
        ('STA_INDEX', True, '>i2', (3,),   _u.is_strictpos, None, None),
        ('U1COORD',   True, '>f8', (),     None,            None, "m"), 
        ('V1COORD',   True, '>f8', (),     None,            None, "m"),
        ('U2COORD',   True, '>f8', (),     None,            None, "m"), 
        ('V2COORD',   True, '>f8', (),     None,            None, "m"),
        ('T3AMP',     True, '>f8', (_NW,), None,            None, None), 
        ('T3AMPERR',  True, '>f8', (_NW,), None,            None, None),
        ('T3PHI',     True, '>f8', (_NW,), None,            None, "deg"), 
        ('T3PHIERR',  True, '>f8', (_NW,), None,            None, "deg"),
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
        ('CORRINDX_T3AMP', False, '>i4', (), _u.is_strictpos, None, None), 
        ('CORRINDX_T3PHI', False, '>i4', (), _u.is_strictpos, None, None),
    ]

