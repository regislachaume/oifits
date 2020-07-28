from .data import _DataHDU11, _DataHDU22
from .t2 import _T2HDU
from .. import utils as _u

import numpy as np

class _VisHDU(_T2HDU):
    _EXTNAME = 'OI_VIS'    
    _COLUMNS = [
        ('VISAMP',    True, '>f8', ('NWAVE',), None, None, None), 
        ('VISAMPERR', True, '>f8', ('NWAVE',), None, None, None),
        ('VISPHI',    True, '>f8', ('NWAVE',), None, None, None), 
        ('VISPHIERR', True, '>f8', ('NWAVE',), None, None, None),
    ]


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
        ('CORRINDX_VISAMP',  True, '>i4', (),           _spos, None, None), 
        ('CORRINDEX_VISPHI', True, '>i4', (),           _spos, None, None),
        ('RVIS',             True, '>f8', ('NWAVE',),   None,  None, None), 
        ('RVISERR',          True, '>f8', ('NWAVE',),   None,  None, None),
        ('IVIS',             True, '>f8', ('NWAVE',),   None,  None, None), 
        ('IVISERR',          True, '>f8', ('NWAVE',),   None,  None, None),
        ('CORRINDX_RVIS',    True, '>i4', (),           _spos, None, None), 
        ('CORRINDEX_IVIS',   True, '>i4', (),           _spos, None, None),
        ('VISREFMAP',        True, '|b1', ('NWAVE',)*2, None,  None, None),
    ]
    
    _CARDS = [
        ('AMPTYP',   False, _is_amptyp2,  None),
        ('PHITYP',   False, _is_phityp2,  None),
        ('AMPORDER', False, _u.is_posint, 0),
        ('PHIORDER', False, _u.is_posint, 0),
    ]
    
    def get_reference_channel(self, output_dim='data', flatten=False):

        visrefmap = getattr(self, 'VISREFMAP', 0)
        if visrefmap == 0:
            return super().get_reference_channel('data', False)

        nwave = visrefmap.shape[-1]
        visrefmap = np.array(visrefmap, dtype=object) # allows long int's
        visrefmap = sum(visrefmap[:,:,i] << (i + 1) for i in range(nwave))
        if flatten:
            visrefmap = visrefmap.ravel()
 
        return visrefmap
    
    def get_obs_type(self, name, output_dim='data', flatten=False):
        
        if name == 'VISAMP':
            typ = self.header.get('AMPTYP', 'N/A')
        elif name == 'VISPHI':
            typ = self.header.get('PHITYP', 'N/A')
        else:
            typ = 'correlated flux'
        return self._resize_data(typ, output_dim, flatten)
