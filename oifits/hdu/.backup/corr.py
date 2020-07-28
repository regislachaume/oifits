from .base import _OIExtHDU, _OIExtHDU1
from .. import utils as _u

class _CorrBase(_OiHDU):
    
    def get_corrname(self, output_dim='none', flatten=False, default=None):
        corrname = self.header.get('CORRNAME', default)
        return self._resize_data(corrname, output_dim, flatten)

    def get_corrHDU(self):
        corrname = self.get_corrname()
        if not corrname:
            return None
        return self.container.get_corrHDU(corrname)

class _MayHaveCorrHDU(_CorrBase):
    _CARDS = [('CORRNAME', False, _u.is_nonempty, None)]

class _MustHaveCorrHDU(_CorrBase):
    _CARDS = [('CORRNAME', True,  _u.is_nonempty, None)]

class _CorrHDU(_MustHaveCorrHDU):
    _EXTNAME = 'OI_CORR'

class CorrHDU1(
        _CorrHDU,
        _OiHDU1,
        _BinTableHDU,
    ):
    _COLUMNS = [
        ('IINDX', True, '<i4', (), _u.is_strictpos, None, None),
        ('JINDX', True, '<i4', (), _u.is_strictpos, None, None), 
        ('CORR',  True, '<f8', (), None,            None, None)
    ]
