from .table import _OITableHDU
from .. import utils as _u

class _T2HDU(_OITableHDU):
    
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
    
