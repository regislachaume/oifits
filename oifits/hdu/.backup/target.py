import numpy as np

from .base import _OIExtHDU, _OiHDU1, _OiHDU2
from .. import utils as _u

__all__ = ["TargetHDU1", "TargetHDU2"]

milliarcsec = np.deg2rad(1) / 3_600_000

class _MustHaveTargetHDU(_OIExtHDU):

    def _get_target_field(self, name, output_dim='none', flatten=False,
            default=None):

        refhdu = self.get_targetHDU()
        if refhdu is None:
            val = default
        else:
            val = base._xmatch(self, name, refhdu, 'TARGET_ID')
        return self._resize_data(val, output_dim, flatten)

    def get_target(self, output_dim='none', flatten=False, default='N/A'):
        obj = self.container[0].header.get('OBJECT', default)
        if obj != 'MULTI':
            default = obj
        return self._get_target_field('TARGET', output_dim, flatten, default)
    def get_equinox(self, output_dim='none', flatten=False):
        return self._get_target_field('EQUINOX', output_dim, flatten)
    def get_ra(self, output_dim='none', flatten=False):
        return self._get_target_field('RAEP0', output_dim, flatten)
    def get_dec(self, output_dim='none', flatten=False):
        return self._get_target_field('DECEP0', output_dim, flatten)
    def get_parallax(self, output_dim='none', flatten=False):
        return self._get_target_field('PARALLAX', output_dim, flatten)
    def get_pmra(self, output_dim='none', flatten=False):
        return self._get_target_field('PMRA', output_dim, flatten)
    def get_pmdec(self, output_dim='none', flatten=False):
        return self._get_target_field('PMDEC', output_dim, flatten)
    def get_rv(self, output_dim='none', flatten=False):
        return self._get_target_field('SYSVEL', output_dim, flatten)
    def get_spectype(self, output_dim='none', flatten=False):
        return self._get_target_field('SPECTYP', output_dim, flatten)
    def get_category(self, output_dim='none', flatten=False):
        return self._get_target_field('CATEGORY', output_dim, flatten)
    
    def get_targetHDU(self):
        return self.container.get_targetHDU()

    def _verify(self, option='warn'):

        errors = super()._verify(option)

        # Verify Target ID is correct (> 0) and referenced

        val = np.unique(self.TARGET_ID)
        if not all(val >= 1):
            err_text = "'TARGET_ID' should be ≥ 1"
            err = self.run_option(option, err_text, fixable=False)
            errors.append(err)

        t = self.get_targetHDU()
        if t is not self:
            for v in val:
                if v not in t.TARGET_ID:
                    err_text = f"'TARGET_ID' not refered in TargetHDU: {v}"
                    err = self.run_option(option, err_text, fixable=False)
                    errors.append(err)

        return errors


def _is_veltyp(s):
    return s in ['LSR', 'HELIOCEN', 'BARYCENTR', 'TOPOCENT']

def _is_veldef(s):
    return s in ['OPTICAL', 'RADIO']

class _TargetHDU(_MustHaveTargetHDU):
    _EXTNAME = 'OI_TARGET' 
    _COLUMNS = [
        ('TARGET_ID',  True, '>i2',  (), _u.is_strictpos, None, None),
        ('RAEP0',      True, '>f8',  (), None,            None, "deg"), 
        ('DECEP0',     True, '>f8',  (), None,            None, "deg"), 
        ('EQUINOX',    True, '>f4',  (), None,            None, "yr"),
        ('RA_ERR',     True, '>f8',  (), None,            0.,   "deg"),  
        ('DEC_ERR',    True, '>f8',  (), None,            0.,   "deg"),
        ('SYSVEL',     True, '>f8',  (), None,            None, "m/s"), 
        ('VELTYP',     True, '|S8',  (), _is_veltyp,      None, None), 
        ('VELDEF',     True, '|S8',  (), _is_veldef,      None, None),
        ('PMRA',       True, '>f8',  (), None,            0.,   "deg/yr"), 
        ('PMDEC',      True, '>f8',  (), None,            0.,   "deg/yr"),
        ('PMRA_ERR',   True, '>f8',  (), None,            0.,   "deg/yr"), 
        ('PMDEC_ERR',  True, '>f8',  (), None,            0.,   "deg/yr"),
        ('PARALLAX',   True, '>f4',  (), None,            None, "deg"), 
        ('PARA_ERR',   True, '>f4',  (), None,            0.,   "deg"),
    ]
    
    def _verify(self, option='warn'):

        errors = super()._verify(option)

        val = getattr(self, 'CATEGORY', None)
        if val is not None and not all((val == 'SCI') * (val == 'CAL')):
            err_text = "'CATEGORY' should be either SCI or CAL"
            err = self.run_option(option, err_text, fixable=False)
            errors.append(err)

        return errors
    
    def _argwhere_same_target(self, t2, max_distance=10 * milliarcsec):
        indices = np.logical_and(
                      self.TARGET == t2['TARGET'],
                      self.EQUINOX == t2['EQUINOX'],
                      abs(self.RAEP0 - t2['RAEP0']) < max_distance,
                      abs(self.DECEP0 - t2['DECEP0']) < max_distance
                  )
        return np.argwhere(indices)
                
    def _merge(self, hdu2, max_distance=10 * milliarcsec):

        # re-index first HDU
        old_id1 = self.TARGET_ID
        new_id1 = np.arange(1, 1 + len(old_id1)) 
        
        # append second HDU, but only lines that are different
        i = 1 + len(old_id1)
        old_id2 = hdu2.TARGET_ID
        new_id2 = np.zeros_like(old_id2)
        kept_lines = []
        for j, t2 in enumerate(hdu2.data):
            where = hdu1._argwhere_same_target(t2, max_distance=max_distance)
            if len(where):
                new_id2[j] = old_id1[where[0,0]] 
            else:
                new_id2[j] = i
                i += 1
                kept_lines.append(j)
        
        # merged table
        merged = self + hdu2.data[kept_lines] 
        merged.TARGET_ID = np.hstack(new_id1, new_id2[kept_lines])
        
        # maps old to new indices for each table
        map1 = {o: n for o, n in zip(old_id1, new_id1)}
        map2 = {o: n for o, n in zip(old_id2, new_id2)}
        
        return merged, map1, map2
 
def _is_category(s):
    return s in ['SCI', 'CAL']

class TargetHDU1(
        _TargetHDU,
        _OIExtHDU1,
      ):
    _COLUMNS = [
        ('TARGET',  True, '|S16', (), _u.is_nonempty, None, None),
        ('SPECTYP', True, '|S16', (), None,           None, None),
    ]

class TargetHDU2(
        _TargetHDU,
        _OIExtHDU2,
      ):
    _COLUMNS = [
        ('TARGET',   True,  '|S32', (), _u.is_nonempty, None,  None),
        ('SPECTYP',  True,  '|S32', (), None,           None,  None),
        ('CATEGORY', False, '|S3',  (), _is_category,   'SCI', None),
    ]
