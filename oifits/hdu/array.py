
from .table import _OITableHDU, _OITableHDU11, _OITableHDU22
from .. import utils as _u

__all__ = ["ArrayHDU1", "ArrayHDU2"]

import numpy as _np

class _MayHaveArrayHDUBase(_OITableHDU):

    def _get_array_field(self, name, shape='none', flatten=False,
            concatenate=False, default=None):

        refhdu = self.get_arrayHDU()

        if refhdu is None:
            if default is None:
                return None
            val = default
        else:
            val = self._xmatch(name, refhdu, 'STA_INDEX', concatenate)

        return self._resize_data(val, shape, flatten)

    def get_sta_name(self, shape='none', flatten=False):
        return self._get_array_field('STA_NAME', shape, flatten)

    def get_sta_xyz(self, shape='none', flatten=False):
        return self._get_array_field('STAXYZ', shape, flatten)

    def get_sta_config(self, shape='none', flatten=False, default=None):
        return self._get_array_field('STA_NAME', shape, flatten,
                    concatenate=True, default=default)

    def get_tel_name(self, shape='none', flatten=False):
        return self._get_array_field('TEL_NAME', shape, flatten)

    def get_tel_config(self, shape='none', flatten=False, default=None):
        return self._get_array_field('TEL_NAME', shape, flatten,
                    concatenate=True, default=default)

    def get_arrname(self, shape='none', flatten=False, default=None):
        arrname = self.header.get('ARRNAME', default)
        if not arrname:
            return None
        arrname = self._resize_data(arrname, shape, flatten)
        return arrname

    def get_arrayHDU(self):
        return self._container.get_arrayHDU(arrname=self.get_arrname())

    def _verify(self, option='warn'):

        errors = super()._verify(option)
      
        # If no ARRNAME and STA_INDEX, fine (optional), but if
        # only one of them is present, there's a problem.
        sta_index = getattr(self, 'STA_INDEX', None)
        arrname = self.header.get('ARRNAME', None) 
        
        if not arrname:
            if sta_index:
                err_text = "STA_INDEX column present with no ARRNAME keyword"
                err = self.run_option(option, err_text, fixable=False)
                errors.appen(err)
            return errors
         
        if sta_index is None:
            err_text = "ARRNAME is given but not STA_INDEX column"
            err = self.run_option(option, err_text, fixable=False)
            errors.append(err)
            return errors

        # check it has the corresponding arrayHDU
        refhdu = self.get_arrayHDU() 
        if refhdu is None:
            err_text = "ArrayHDU with ARRNAME={arrname} not found"
            err = self.run_option(option, err_text, fixable=False)
            errors.append(err)
            return errors

        # We don't report on i <= 0 unreferenced indices: they are
        # already spotted as not valid by column verification 
        ref_index = refhdu.STA_INDEX
        sta_index = _np.unique(sta_index)
        missing = [i for i in sta_index if i > 0 and i not in ref_index]
        if len(missing):
            m = ', '.join(missing)
            err_text = f"'STA_INDEX' not referenced in ArrayHDU: {m}"
            err = self.run_option(option, err_text, fixable=False)
            errors.append(err)

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
        ('TEL_NAME',  True, '<U8', (),   None,            None, None), 
        ('STA_NAME',  True, '<U8', (),   None,            None, None), 
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

    def _verify(self, option='warn'):

        errors = super()._verify(option)
        
        sta_index = self.STA_INDEX
        if len(_np.unique(sta_index)) == len(sta_index):
            return errors

        err_text = f"Repeated STA_INDEX in {type(self).__name__}"
        err = self.run_option(option, err_text, fixable=False)
        errors.append(err)

        return errors

    def _merge(self, other, max_distance=0.05):

        # how we decide two stations are the same
        id_key = 'STA_INDEX'
        eq_keys = ['STA_NAME', 'TEL_NAME'] # must be equal 
        dist_keys = ['STAXYZ'] # may differ by 5 cm
        
        return super()._merge(other, id_key, eq_keys, dist_keys,
                                max_distance=max_distance)


class ArrayHDU1(
        _ArrayHDU,
        _OITableHDU11, # OFITS1, rev. 1
      ):
    _CARDS = [('FRAME', True, _u._is_frame1, 'GEOCENTRIC')]


class ArrayHDU2(
        _ArrayHDU,
        _OITableHDU22, # OIFITS2, rev. 2
      ):
    _CARDS = [('FRAME', True, _u._is_frame2, 'GEOCENTRIC')]
    _COLUMNS = [
        ('FOV',     False, '>f8', (), _u.is_strictpos, None, "arcsec"), 
        ('FOVTYPE', False, '<U6', (), _u._is_fovtype2, None, None)
    ]
