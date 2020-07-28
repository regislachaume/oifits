from numpy import array

from .base import _OIExtHDU, _OIExtHDU1, _OIExtHDU2
from .. import utils as _u

__all__ = ["ArrayHDU1", "ArrayHDU2"]

class _MayHaveArrayHDUBase(_OIExtHDU):

    def _get_array_field(self, name, output_dim='none', flatten=False,
            concatenate=False, default=None):

        refhdu = self.get_arrayHDU()

        if refhdu is None:
            if default is None:
                return None
            val = default
        else:
            val = base._xmatch(self, name, refhdu, 'STA_INDEX', concatenate)

        return self._resize_data(val, output_dim, flatten)

    def get_sta_name(self, output_dim='none', flatten=False):
        return self._get_array_field('STA_NAME', output_dim, flatten)

    def get_sta_xyz(self, output_dim='none', flatten=False):
        return self._get_array_field('STAXYZ', output_dim, flatten)

    def get_sta_config(self, output_dim='none', flatten=False, default=None):
        return self._get_array_field('STA_NAME', output_dim, flatten,
                    concatenate=True, default=default)

    def get_tel_name(self, output_dim='none', flatten=False):
        return self._get_array_field('TEL_NAME', output_dim, flatten)

    def get_tel_config(self, output_dim='none', flatten=False, default=None):
        return self._get_array_field('TEL_NAME', output_dim, flatten,
                    concatenate=True, default=default)

    def get_arrname(self, output_dim='none', flatten=False, default=None):
        arrname = self.header.get('ARRNAME', default)
        if not arrname:
            return None
        arrname = self._resize_data(arrname, output_dim, flatten)
        return arrname

    def get_arrayHDU(self):
        return self.container.get_arrayHDU(arrname=self.get_arrname())

    def _verify(self, option='warn'):

        errors = super()._verify(option)
        
        # only non ArrayHDU needs to be checked for cross references
        arrname = self.header.get('ARRNAME', None) 
        if isinstance(self, _ArrayHDU) or not arrname:
            return errors

        # check it has the corresponding arrayHDU
        arrayHDU = self.get_arrayHDU() 
        ok = arrayHDU is not None
        self.req_cards('ARRNAME', None, lambda s: ok, None, option, errors)
        if not ok:
            return errors
        
        # check the cross references
        indices = np.unique(self.STA_INDEX)
        ref_indices = arrayHDU.STA_INDEX
        ok = all(i in ref_indices for i in indices)

        # check it 
        return errors

class _MayHaveArrayHDU(_MayHaveArrayHDUBase):
    _CARDS = [('ARRNAME', False, _u.is_nonempty, None)]

class _MustHaveArrayHDU(_MayHaveArrayHDUBase):
    _CARDS = [('ARRNAME', True, _u.is_nonempty, None)] 

class _ArrayHDU(_MustHaveArrayHDU):
    
    _EXTNAME = 'OI_ARRAY'
    _CARDS = [
        ('ARRAYX', True, _u.is_num, None),
        ('ARRAYY', True, _u.is_num, None),
        ('ARRAYZ', True, _u.is_num, None),
    ]
    _COLUMNS = [
        ('TEL_NAME',  True, '|S8', (),   None,            None, None), 
        ('STA_NAME',  True, '|S8', (),   None,            None, None), 
        ('STA_INDEX', True, '>i2', (),   _u.is_strictpos, None, None),
        ('DIAMETER',  True, '>f4', (),   _u.is_strictpos, None, "m"),
        ('STAXYZ',    True, '>f8', (3,), None,            None, "m"),
    ]
     
    def _get_ins(self):
        name = self.arrname()
        m = re.match('[A-Za-z_]{1,20}', name)
        if m:
            return m.group()
        return name[0:20] 

def _is_frame1(s):
    return s == 'GEOCENTRIC'

class ArrayHDU1(
        _ArrayHDU,
        _OIExtHDU1,
      ):
    _CARDS = [('FRAME', True, _is_frame1, 'GEOCENTRIC')]

def _is_frame2(s):
    return s in ['GEOCENTRIC', 'SKY']
def _is_fovtype2(s):
    return s in ['FVWM', 'RADIUS']

class ArrayHDU2(
        _ArrayHDU,
        _OIExtHDU2,
      ):
    _CARDS = [('FRAME', True, _is_frame2, 'GEOCENTRIC')]
    _COLUMNS = [
        ('FOV',     False, '>f8', (), _u.is_strictpos, None, "arcsec"), 
        ('FOVTYPE', False, '|S6', (), _is_fovtype2,    None, None)
    ]
