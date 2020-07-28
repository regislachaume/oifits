from .data import _DataHDU11, _DataHDU22
from .t2 import _T2HDU
from .. import utils as _u

class _Vis2HDU(_T2HDU):
    _EXTNAME = 'OI_VIS2'
    _COLUMNS = [
        ('VIS2DATA', True, '>f8', ('NWAVE',), None, None, None), 
        ('VIS2ERR',  True, '>f8', ('NWAVE',), None, None, None),
    ]
    
    def get_obs_type(self, name, output_dim='data', flatten=False):
        
        return self._resize_data('absolute', output_dim, flatten)

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
        ('CORRINDX_VIS2DATA', True, '>i4', (), _u.is_strictpos, None, None),
    ]
