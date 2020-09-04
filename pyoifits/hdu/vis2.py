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
    def from_data(cls, insname, *, version=2, arrname=None, corrname=None,
        target_id=None, mjd=None, int_time=_np.nan, time=0., ucoord=0.,
        vcoord=0., vis2data=None, vis2err=None, flag=False,
        fits_keywords={}, **columns):

        assert vis2data is not None, 'vis2data must be specified'
        assert vis2err is not None, 'vis2err must be specified'

        if vis2data is not None:
            shape = _np.shape(vis2data)
            if not _np.shape(ucoord):
                ucoord = _np.full(shape, ucoord)
            if not _np.shape(vcoord):
                vcoord = _np.full(shape, vcoord)
            if not _np.shape(flag):
                flag = _np.full(shape, flag)
        
        return super().from_data(insname, version=2, arrname=arrname, 
            corrname=corrname, target_id=target_id, mjd=mjd, int_time=int_time, 
            time=time, ucoord=ucoord, vcoord=vcoord, vis2data=vis2data, 
            vis2err=vis2err, flag=flag, fits_keywords=fits_keywords, **columns)

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
