from .data import _DataHDU11, _DataHDU22
from .t2 import _T2HDU
from .. import utils as _u
from .wavelength import _NW

import numpy as _np

class _Vis2HDU(_T2HDU):
    _EXTNAME = 'OI_VIS2'
    _COLUMNS = [
        ('VIS2DATA', True, '>f8', (_NW,), None, None, None,
            'squared visibility amplitude'), 
        ('VIS2ERR',  True, '>f8', (_NW,), None, None, None,
            'uncertaintly on squared visibility amplitude'),
    ]
    
    def get_obs_type(self, name, shape='data', flatten=False):
        
        return self._resize_data('absolute', shape, flatten)

    @classmethod
    def from_data(cls, *, vis2data, fits_keywords={}, **columns):

        columns = dict(vis2data=vis2data, **columns)

        shape = self._get_columns_shape(**columns)
        _u.store_default(columns, 'vis2err', 0., shape)

        return super().from_data(insname, arrname=arrname, 
            corrname=corrname, version=2, date=date,
            fits_keywords=fits_keywords, **columns)

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
        ('CORRINDX_VIS2DATA', False, '>i4', (), _u.is_strictpos, None, None,
            'index of 1st visib. in matching OI_CORR matrix'),
    ]

new_vis2_hdu = _Vis2HDU.from_data
