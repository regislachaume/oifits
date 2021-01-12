from .array import _MustHaveArrayHDU
from .table import _OITableHDU, _OITableHDU21
from .wavelength import _NW
from .. import utils as _u

import numpy as _np

class _InspolHDUBase(_OITableHDU):

    def get_inspolHDU(self):

        if isinstance(self, InspolHDU1):
            return self

        arrname = self.header.get('ARRNAME')
        def same_arrname(h): return h.get('ARRNAME', '') == arrname
        container = self.get_container()
        
        hdus = container.get_HDUs(_InspolHDU, filter=same_arrname)
        if not hdus:
            return None
        return hdus[0]

    def get_jones_matrix(self, shape='data', flatten=False):

        J = _np.zeros((*self.data_shape(), 2, 2), dtype=complex)

        inspol = self.get_inspolHDU()

        if inspol is not None:

            insname = self.get_insname()
            poldata = inspol.data

            for i, row in enumerate(self.data):
                match = ((poldata['INSNAME'] == insname)
                        * (poldata['TARGETID'] == row['TARGETID'])
                        * (poldata['MJD'] <= row['MJD'])
                        * (row['MDJ'] <= poldata['MJD-END']))
                if any(match):
                    j = poldata[match][0]
                    J[i,:,0,0] = j['JXX']
                    J[i,:,0,1] = j['JXY']
                    J[i,:,1,0] = j['JYX']
                    J[i,:,1,1] = j['JYY']

        if flatten:
            J = J.reshape((-1, 2, 2))

        return J

class _MayHaveInspolHDU(
        _InspolHDUBase,
        # _MayHaveArrayHDU, -> would give incorrect card
      ):
    pass

class _MustHaveInspolHDU(
        _InspolHDUBase,
        _MustHaveArrayHDU,
      ):
    pass 
    
class _InspolHDU(
        _MustHaveInspolHDU,
      ):
    
    _EXTNAME = 'OI_INSPOL'
    
    # exactly similar to _DataHDU.  Rethink class hierarchy?
    def _trim_helper(self, target_filer=None,    
        wave_filter=lambda wave: True, insname_filter=lambda ins: True,
        keep_ns_columns=False):

        insname = self.get_insname()
        whdu = self.get_wavelength_HDU(insname=insname[0])
        wkeep = _np.vectorize(wave_filter)(whdu.get_wave())
        tkeep = _np.vectorize(target_filter)(self.get_target())
        ikeep = _np.vectorize(insname_filter)(insname)

        columns = {}
        spectral_colnames = self._get_spec_colnames()
        standard_colnames = self._get_oi_colnames()
        for column in self.columns:
            colname = column.name
            if not keep_ns_columns and colname not in standard_colnames:
                continue
            data = self.data[colname]
            if colname in spectral_colnames:
                data = data[:,wkeep]
            columns[colname] = data[tkeep & ikeep]

        phdu = self._from_data(fits_keywords=self.header, **columns)

        return phdu
    
class InspolHDU1(
        _InspolHDU, 
        _OITableHDU21, # OIFITS2, table rev. 1
      ):
    """

First revision of the OI_INSPOL binary table, OIFITS v. 2

    """
    _CARDS = [
        ('NPOL',   True, _u.is_strictpos, None, 
            'number of polarisations'),
        ('ORIENT', True, _u.is_nonempty,  None, 
            'orientation of the Jones matrix'),
        ('MODEL',  True, _u.is_nonempty,  None, 
            'method used to determine the Jones matrix')
    ]
    _COLUMNS = [
        ('TARGET_ID',  True, '1I',  (),     _u.is_strictpos, None, None,
            'target ID in matching OI_TARGET table'),
        ('INSNAME',    True, '32A', (),     None,            None, None,
            'name of matching OI_WAVELENGTH table'), 
        ('MJD_OBS',    True, '1D',  (),     None,            None, "d",
            'modified Julian day at start of observation'), 
        ('MJD_END',    True, '1D',  (),     None,            None, "d",
            'modified Julian day at en of observation'),
        ('JXX',        True, 'M',   (_NW,), None,            None, None,
            'Jones matrix element J_xx'), 
        ('JYY',        True, 'M',   (_NW,), None,            None, None,
            'Jones matrix element J_xy'),
        ('JXY',        True, 'M',   (_NW,), None,            None, None,
            'Jones matrix element J_yx'), 
        ('JYX',        True, 'M',   (_NW,), None,            None, None,
            'Jones matrix element J_yy'),
        ('STA_INDEX',  True, '1I',  (),     None,            None, None,
            'station index in matching OI_ARRAY table'),
    ]

    def __mod__(self, other):
        return False

    def get_insname(self, shape='data', flatten=False):

        insname = self.data['INSNAME']
        return self._resize_data(insname, shape, flatten)

    def get_jones_matrix(self, shape='data', flatten=False):

        poldata = self.data

        nrows, nwaves = poldata['JXX'].shape
        J = _np.empty((nrows, nwaves, 2, 2), dtype=complex)

        J[...,0,0] = poldata['JXX']
        J[...,0,1] = poldata['JXY']
        J[...,1,0] = poldata['JYX']
        J[...,1,1] = poldata['JYY']

        if flatten:
            J = J.reshape((-1, 2, 2))

        return J

