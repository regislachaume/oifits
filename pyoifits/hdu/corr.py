"""
Implementation of the OI_CORR binary table extension containing corralation
between interferometric data in OI_VIS, OI_VIS2, OI_T3, and OI_FLUX binary
table extensions
"""

from .table import _OITableHDU, _OITableHDU21
from .referenced import _Referenced
from .. import utils as _u
from scipy import sparse as _sparse

import numpy as _np

class _CorrHDUBase(_OITableHDU):
    
    def get_corrname(self, shape='none', flatten=False, default=None):
        corrname = self.header.get('CORRNAME', default)
        return self._resize_data(corrname, shape, flatten)

    def get_corrHDU(self):
        corrname = self.get_corrname()
        if not corrname:
            return None
        return self._container.get_corrHDU(corrname)

class _MayHaveCorrHDU(_CorrHDUBase):
    _CARDS = [('CORRNAME', False, _u.is_nonempty, None, 
        'correlation name for cross-reference')]

class _MustHaveCorrHDU(_CorrHDUBase):
    _CARDS = [('CORRNAME', True,  _u.is_nonempty, None, 
        'correlation name for cross-reference')]

class _CorrHDU(_MustHaveCorrHDU,_Referenced):
    
    _EXTNAME = 'OI_CORR'
    _REFERENCE_KEY = 'CORRNAME'

    _CARDS = [
        ('NDATA', True,  _u.is_strictpos, None, 
            'size of correlation matrix (NDATAxNDATA)'),
    ] 
    _COLUMNS = [
        ('IINDX', True, '1J', (), _u.is_strictpos, None, None,
            '1st index (<= NDATA)'),
        ('JINDX', True, '1J', (), _u.is_strictpos, None, None,
            '2nd index (<= NDATA)'), 
        ('CORR',  True, '1D', (), None,            None, None,
            'correlation C[IINDX,JINDX]')
    ]

    def _trim_helper(self, *, keep_ns_columns=False, 
            wave_filter=None, target_filter=None, insname_filter=None):

        raise NotImplementedError('cannot trim correlated data...')
    
    def _verify(self, option='warn'):

        errors = super()._verify(option)

        corr_index = self.CORR_INDEX
        if len(_np.unique(corr_index)) == len(coor_index):
            return errors

        err_text = f"Repeated CORR_INDEX in {type(self).__name__}"
        err = self.run_option(option, err_text, fixable=False)
        errors.append(err)

        return errors

    def get_corrname(self):
        return self.header['CORRNAME']

    def __mod__(self, other):

        return False

    def is_referred_to_by(self, other):
        return (not isinstance(other, _CorrHDU) and
                isinstance(other, _CorrHDUBase) and
                self.get_corrname() == other.get_corrname())

    @classmethod
    def from_data(cls, *, corrname, corrmatrix, fits_keywords={}, **columns):

        fits_keywords = dict(corrname=corrname, **fits_keywords)

        # find non-zero elements in the lower triangle, note that
        # indices start at one in OIFITS        
        iindx, jindx, corr = _sparse.find(corrmatrix)
        keep = iindx < jindx
        corr = corr[keep]
        iindx = iindx[keep]
        jindx = jindx[keep]
        iindx += 1
        jindx += 1

        columns = dict(iindx=iindx, jindx=jindx, corr=corr, **columns)        

        super().from_data(fits_keywords=fits_keywords, **columns)
    
class CorrHDU1(
        _CorrHDU,
        _OITableHDU21, # OIFITS2, table rev. 1
    ):
    """

First revision of the OI_CORR binary table, OIFITS v. 2

    """
    pass


