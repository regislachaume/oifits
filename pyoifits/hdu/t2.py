from .data import _DataHDU
from .. import utils as _u

class _T2HDU(_DataHDU):
    
    _COLUMNS = [
        ('STA_INDEX', True, '>i2', (2,), _u.is_strictpos, None, None,
            'station indices in matching OI_ARRAY table'),
        ('UCOORD',    True, '>f8', (),   None,            None, "m",
            'u coordinate'),
        ('VCOORD',    True, '>f8', (),   None,            None, "m",
            'v coordinate'),
    ]
    
    def __getattr__(self, name):

        if name == 'U1COORD':
            return self.data['UCOORD']
        if name == 'V1COORD':
            return self.data['VCOORD']
        
        return super().__getattr__(name)
    
