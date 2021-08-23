import numpy as np

from .table import _OITableHDU, _OITableHDU21
from .referenced import _Referenced
from .. import utils as u

__all__ = ["CorrHDU1"]

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

    _CARDS = [('CORRNAME', False, u.is_nonempty, None, 
        'correlation name for cross-reference')]
    
    def get_corrindx(self, obsname=None, shape='none', flatten=False):

        if obsname is None:
            colnames = self.columns.names
            obsnames = self.get_observable_names()
            obsnames = [n for n in obsnames if n in colnames]
            index = [self.get_corrindx(obsname, flatten=flatten, shape=shape)
                            for obsname in obsnames]
            return index
        
        corrindex_name = obsname + '_CORRINDEX'
        if corrindex_name not in self.columns.names:
            index = np.zeros_like(self.data[obsname], dtype=int)
            index = np.ma.masked_equal(index, 0)
        else:
            index = self.data[corrindex_name]
            relindex = np.arange(self.get_nwaves())
            index = index[:,None] + relindex
            index = np.ma.masked_array(index, mask=self.FLAG)

        if flatten:
            index = index.ravel()

        return index

    def _verify(self, option='warn'):
        
        errors = super()._verify(option)
   
        # check for unicity of indices
 
        index = np.ravel([i[~i.mask] for i in self.get_corrindx(flatten=True)])
        
        if len(np.unique(index)) < len(index):
            err_text = f"repeated CORRINDX in {self.header['EXTNAME']}"
            err = self.run_option(option, err_text, fixable=False)
            errors.append(err)

        return errors

class _CorrHDU(_CorrHDUBase,_Referenced):
    
    _EXTNAME = 'OI_CORR'
    _REFERENCE_KEY = 'CORRNAME'

    _CARDS = [
        ('CORRNAME', True,  u.is_nonempty, None, 
            'correlation name for cross-reference'),
        ('NDATA', True,  u.is_strictpos, None, 
            'size of correlation matrix (NDATAxNDATA)'),
    ] 
    _COLUMNS = [
        ('IINDX', True, '1J', (), u.is_strictpos, None, None,
            '1st index'),
        ('JINDX', True, '1J', (), u.is_strictpos, None, None,
            '2nd index'), 
        ('CORR',  True, '1D', (), None,            None, None,
            'correlation C[IINDX,JINDX]')
    ]

    def _trim_helper(self, *, keep_ns_columns=False, 
            wave_filter=None, target_filter=None, insname_filter=None):

        raise NotImplementedError('cannot trim correlated data...')
    
    def _verify(self, option='warn'):

        errors = super()._verify(option)

        # order IINDX JINDX

        inverted = self.IINDX > self.JINDX
        if any(inverted):
            err_text = "IINDX > JINDX in OI_CORR."
            fix_text = "Order swapped"
            def fix(h=self):
                tmp = self.IINDX[inverted]
                self.IINDX[inverted] = self.JINDX[inverted]
                self.JINDX[inverted] = tmp
            err = self.run_option(option, err_text, fix_text, fix)

        # ignore IINDX = JINDX

        not_equal = self.IINDX != self.JINDX
        if not all(not_equal):
            err_text = 'IINDX = JINDX in OI_CORR.'
            fix_text = 'Removed'
            def fix(h=self):
                h.data = h.data[not_equal]
                h.header['NAXIS2'] = len(h.data)
            err = self.run_option(option, err_text, fix_text, fix) 

        # check for redundancies
        n = len(self.IINDX)
        n_unique = np.unique([self.IINDX, self.JINDX], axis=1).shape[1]

        if n_unique < n:
            err_text = f"Repeated CORR_INDEX pair in OI_CORR."
            err = self.run_option(option, err_text, fixable=False)
            errors.append(err)

        # check for ambiguous references from data HDUs
        hdus = self.get_referrers()
        index = np.ma.ravel([h.get_corrindx(flatten=True) for h in hdus])
        index = index[~index.mask]
        if len(np.unique(index)) < len(index):
            err_text = f"Different data refer to same OI_CORR matrix element"
            err = self.run_option(option, err_text, fixable=False)
            errors.append(err)
        
        # check for matrix elements not referred to
        unref = [i not in index for i in np.unique([self.IINDX, self.JINDX])]
        if len(unref):
            err_text = f"OI_CORR matrix element(s) not referred to by any data"
            fix_text = f"Removed element(s)"
            keep = ~np.in1d(self.IINDX, unref) & ~np.in1d(self.JINDX, unref)
            def fix(h=self):
                h.data = h.data[keep]
            err = self.run_option(option, err_text, fix_text, fix)
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
    
        from scipy import sparse

        fits_keywords = dict(corrname=corrname, **fits_keywords)

        # find non-zero elements in the lower triangle, note that
        # indices start at one in OIFITS        
        iindx, jindx, corr = sparse.find(corrmatrix)
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


