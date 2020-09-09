from .data import _DataHDU11, _DataHDU22
from .t2 import _T2HDU
from .wavelength import _NW
from .. import utils as _u

import numpy as np

class _VisHDU(_T2HDU):
    _EXTNAME = 'OI_VIS'    
    _COLUMNS = [
        ('VISAMP',    True, '>f8', (_NW,), None, None, None,
            'visibility amplitude'), 
        ('VISAMPERR', True, '>f8', (_NW,), None, None, None,
            'uncertainty on visibility amplitude'),
        ('VISPHI',    True, '>f8', (_NW,), None, None, None,
            'phase'), 
        ('VISPHIERR', True, '>f8', (_NW,), None, None, None,
            'uncertainty on phase'),
    ]
    
    @classmethod
    def from_data(cls, *, visamp, visphi, fits_keywords={}, **columns):

        columns = dict(visamp=visamp, visphi=visphi,  **columns)

        shape = self._get_columns_shape(**columns)
        for name in ['visamp', 'visphi', 'rvis', 'ivis']:
            if name in columns:
                _u.store_default(columns, f"{name}err", 0., shape)

        return super().from_data(fits_keywords=fits_keywords, **columns)

class VisHDU1(
        _VisHDU,
        _DataHDU11,
      ):
    pass

_spos = _u.is_strictpos

def _is_amptyp2(s):
    return s in ['absolute', 'differential', 'correlated flux']

def _is_phityp2(s):
    return s in ['absolute', 'differential']


class VisHDU2(
        _VisHDU,
        _DataHDU22,
      ):
    _COLUMNS = [
        ('CORRINDX_VISAMP',  False, '>i4', (),         _spos, None, None,
            'index on 1st amp. in matching OI_CORR matrix'), 
        ('CORRINDEX_VISPHI', False, '>i4', (),         _spos, None, None,
            'index on 1st phase in matching OI_CORR matrix'),
        ('RVIS',             False, '>f8', (_NW,),     None,  None, None,
            'real part of correlated flux'), 
        ('RVISERR',          False, '>f8', (_NW,),     None,  None, None,
            'uncertainty on real part of correlated flux'),
        ('IVIS',             False, '>f8', (_NW,),     None,  None, None,
            'imag. part of correlated flux'), 
        ('IVISERR',          False, '>f8', (_NW,),     None,  None, None,
            'uncertainty on imag. part of correlated flux'),
        ('CORRINDX_RVIS',    False, '>i4', (),         _spos, None, None,
            'index of 1st RVIS in matching OI_CORR matrix'), 
        ('CORRINDEX_IVIS',   False, '>i4', (),         _spos, None, None,
            'index of 1st IVIS in matching OI_CORR matrix'),
        ('VISREFMAP',        False, '|b1', (_NW,_NW,), None,  None, None,
            'reference channels for differential quantities'),
    ]
    
    _CARDS = [
        ('AMPTYP',   False, _is_amptyp2,  None, 
            'amplitude is absolute or differential'),
        ('PHITYP',   False, _is_phityp2,  None, 
            'phase is absolute or differential'),
        ('AMPORDER', False, _u.is_posint, 0, 
            'polynomial order in diff. amplitude'),
        ('PHIORDER', False, _u.is_posint, 0, 
            'polynomial order in diff. phase'),
    ]
    
    def get_reference_channel(self, shape='data', flatten=False):

        visrefmap = getattr(self, 'VISREFMAP', 0)
        if visrefmap == 0:
            return super().get_reference_channel('data', False)

        nwave = visrefmap.shape[-1]
        visrefmap = np.array(visrefmap, dtype=object) # allows long int's
        visrefmap = sum(visrefmap[:,:,i] << (i + 1) for i in range(nwave))
        if flatten:
            visrefmap = visrefmap.ravel()
 
        return visrefmap
    
    def get_obs_type(self, name, shape='data', flatten=False):
        
        if name == 'VISAMP':
            typ = self.header.get('AMPTYP', 'N/A')
        elif name == 'VISPHI':
            typ = self.header.get('PHITYP', 'N/A')
        else:
            typ = 'correlated flux'
        return self._resize_data(typ, shape, flatten)

new_vis_hdu = _VisHDU.from_data
