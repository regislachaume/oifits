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

    @classmethod
    def from_data(cls, *, fits_keywords={}, **columns):

        shape = self._get_columns_shape(**columns)
        _u.store_default(columns, 'ucoord', 0., (shape[0],))
        _u.store_default(columns, 'vcoord', 0., (shape[0],)) 
        _u.store_default(columns, 'time', 0., (shape[0],))
        _u.store_default(columns, 'int_time', 0., (shape[0],))
         
        return super().from_data(insname=insname, arrname=arrname, 
            corrname=corrname, version=version, date=date,
            fits_keywords=fits_keywords, columns=columns) 
