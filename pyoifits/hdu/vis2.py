from .data import _DataHDU11, _DataHDU22
from .t2 import _T2HDU
from .. import utils as _u
from .wavelength import _NW

class _Vis2HDU(_T2HDU):
    _EXTNAME = 'OI_VIS2'
    _COLUMNS = [
        ('VIS2DATA', True, '>f8', (_NW,), None, None, None), 
        ('VIS2ERR',  True, '>f8', (_NW,), None, None, None),
    ]
    
    def get_obs_type(self, name, shape='data', flatten=False):
        
        return self._resize_data('absolute', shape, flatten)

class Vis2HDU1(
        _Vis2HDU,
        _DataHDU11,
      ):
    pass

class Vis2HDU2(
        _Vis2HDU,
        _DataHDU22,
      ):
    _COLUMNS = [
        ('CORRINDX_VIS2DATA', False, '>i4', (), _u.is_strictpos, None, None),
    ]
