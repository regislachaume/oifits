from .base import _OIExtHDU
from .. import utils as _u

class _T2HDU(_OIExtHDU):
    
    _COORD_COLUMNS = ['UCOORD', 'VCOORD']
    _COLUMNS = [
        ('STA_INDEX', True, '>i2', (2,), _u.is_strictpos, None, None),
        ('UCOORD',    True, '>f8', (),   None,            None, "m"),
        ('VCOORD',    True, '>f8', (),   None,            None, "m"),
    ]
    
    def __getattr__(self, name):

        if name == 'U1COORD':
            return self.data['UCOORD']
        if name == 'V1COORD':
            return self.data['VCOORD']
        
        return super().__getattr__(name)
    
