from .table import _OITableHDU21
from .target import _MustHaveTargetHDU
from .array import _MayHaveArrayHDU
from .wavelength import _MustHaveWavelengthHDU, _NW
from .data import _DataHDU
from .. import utils as _u

def _is_calstat(s):
    return s in ['U', 'C']
def _is_fovtype(s):
    return s in ['RADIUS', 'FWHM']

class _FluxHDU(
        _DataHDU,
        _MustHaveTargetHDU,
        _MayHaveArrayHDU,
        _MustHaveWavelengthHDU,
      ):
    _EXTNAME = 'OI_FLUX'

class FluxHDU1(
        _FluxHDU,
        _OITableHDU21,
      ):
    
    _CARDS = [
        ('CALSTAT', True,  _is_calstat, None, 
            'Flux is (C)alibrated or (U)ncalibrated'),
        ('FOV',     False, _u.is_pos,   None, 'Field of view (arcsec)'),
        ('FOVTYPE', False, _is_fovtype, None, 'Type of field of view'),
    ]
    
    _COLUMNS = [
        ('FLUXDATA',          True,  '>f8', (_NW,), None, None, None,
            'Flux'),
        ('FLUXERR',           True,  '>f8', (_NW,), None, None, None,
            'Flux uncertainty'),
        ('STA_INDEX',         False, '>i2', (),     None, None, None,
            'Station index in matching OI_ARRAY table'),
        ('CORRINDX_FLUXDATA', False, '>i4', (),     None, None, None,
            'Index of 1st flux in matching OI_CORR matrix'),
    ]
    
    def get_obs_type(self, name, shape='data', flatten=False):
   
        typ = f"{'un' if self.header['CALSTAT'] == 'U' else ''}calibrated flux" 
            
        return self._resize_data(typ, shape, flatten)
    
    def _verify(self, option='warn'):

        errors = []
    
        cols = self.columns
        err = None
        if 'FLUXDATA' not in cols.names and 'FLUX' in cols.names:
            err_text = "Column FLUX should be named FLUXDATA"
            fix_text = "Renamed"
            def fix(h=self): h.rename_columns(FLUX='FLUXDATA')
            err = self.run_option(option, err_text, fix_text, fix)

        errors = super()._verify(option)
        if err:
            errors.append(err)

    @classmethod
    def from_data(cls, *, insname=None, arrname=None, corrname=None,
        version=None, date=None, calstat='U', fov=0., fovtype='N/A',
        fits_keywords={}, **columns):

        shape = self._get_columns_shape(**columns)
        _u.store_default(columns, 'int_time', 0., (shape[0],))
        _u.store_default(columns, 'fluxerr', 0., shape)        

        if calstat == 'U':
            fits_keywords = dict(calstat=calstat, **fits_keywords)
        else:
            fits_keywords = dict(fov=fov, fovtype=fovtype,
                calstat=calstat, **fits_keywords)

        return super().from_data(insname=insname, arrname=arrname,
            corrname=corrname, version=version, date=date,
            fits_keywords={}, **columns):
 
