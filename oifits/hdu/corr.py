from .table import _OITableHDU, _OITableHDU21
from .. import utils as _u

import numpy as _np

class _CorrBase(_OITableHDU):
    
    def get_corrname(self, shape='none', flatten=False, default=None):
        corrname = self.header.get('CORRNAME', default)
        return self._resize_data(corrname, shape, flatten)

    def get_corrHDU(self):
        corrname = self.get_corrname()
        if not corrname:
            return None
        return self._container.get_corrHDU(corrname)

class _MayHaveCorrHDU(_CorrBase):
    _CARDS = [('CORRNAME', False, _u.is_nonempty, None)]

class _MustHaveCorrHDU(_CorrBase):
    _CARDS = [('CORRNAME', True,  _u.is_nonempty, None)]

class _CorrHDU(_MustHaveCorrHDU):
    
    _EXTNAME = 'OI_CORR'

    def _verify(self, option='warn'):

        errors = super()._verify(option)

        corr_index = self.CORR_INDEX
        if len(_np.unique(corr_index)) == len(coor_index):
            return errors

        err_text = f"Repeated CORR_INDEX in {type(self).__name__}"
        err = self.run_option(option, err_text, fixable=False)
        errors.append(err)

        return errors



class CorrHDU1(
        _CorrHDU,
        _OITableHDU21, # OIFITS2, table rev. 1
    ):
    _COLUMNS = [
        ('IINDX', True, '<i4', (), _u.is_strictpos, None, None),
        ('JINDX', True, '<i4', (), _u.is_strictpos, None, None), 
        ('CORR',  True, '<f8', (), None,            None, None)
    ]
