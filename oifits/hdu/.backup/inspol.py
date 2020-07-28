from .array import _MustHaveArrayHDU
from .base import _OIExtHDU, _OIExtHDU1
from .wavelength import _NW
from .. import utils as _u

class _InspolHDUBase(_OIExtHDU):

    def get_inspolHDU(self):

        if isinstance(self, InspolHDU_v1):
            return self

        arrname = self.header.get('ARRNAME')
        def same_arrname(h): return h.get('ARRNAME', '') == arrname
        hdus = self.get_HDUs(InspolHDU, filter=same_arrname)
        if not hdus:
            return None
        return hdus[0]

    def get_jones_matrix(self, output_dim='data', flatten=False):

        J = np.zeros((self.data_shape(), 2, 2), dtype=complex)

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
        _MustHaveArrayHDU
      ):
    pass 
    
class _InspolHDU(_MustHaveInspolHDU):
    _EXTNAME = 'OI_INSPOL'
    
class InspolHDU1(
        _InspolHDU, 
        _OIExtHDU1,
      ):
    _CARDS = [
        ('NPOL',   True, _u.is_strictpos, None),
        ('ORIENT', True, _u.is_nonempty,  None),
        ('MODEL',  True, _u.is_nonempty,  None)
    ]
    _COLUMNS = [
        ('TARGET_ID',  True, '>i2',  (),     _u.is_strictpos, None, None),
        ('INSNAME',    True, '|S32', (),     None,            None, None), 
        ('MJD_OBS',    True, '>f8',  (),     None,            None, "d"), 
        ('MJD_END',    True, '>f8',  (),     None,            None, "d"),
        ('JXX',        True, '>c8',  (_NW,), None,            None, None), 
        ('JYY',        True, '>c8',  (_NW,), None,            None, None),
        ('JXY',        True, '>c8',  (_NW,), None,            None, None), 
        ('JYX',        True, '>c8',  (_NW,), None,            None, None),
        ('STA_INDEX',  True, '>i2',  (),     None,            None, None),
    ]

    def get_insname(self, output_dim='data', flatten=False):

        insname = self.data['INSNAME']
        return self._resize_data(insname, output_dim, flatten)

    def get_jones_matrix(self, output_dim='data', flatten=False):

        J = np.empty((nrows, nwaves, 2, 2), dtype=complex)
        J[...,0,0] = poldata['JXX']
        J[...,0,1] = poldata['JXY']
        J[...,1,0] = poldata['JYX']
        J[...,1,1] = poldata['JYY']

        if flatten:
            J = J.reshape((-1, 2, 2))

        return J

