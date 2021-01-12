"""
Implementation of the OI_FLUX binary table extension containing either
the calibrated fluxes of the targets or the uncalibrated fluxes measured
at the interferometer aperture. 
"""

from .target import _MustHaveTargetHDU
from .array import _MayHaveArrayHDU
from .wavelength import _MustHaveWavelengthHDU, _NW
from .data import _DataHDU, _DataHDU21
from .. import utils as _u


class _FluxHDU(
        _DataHDU,
        _MustHaveTargetHDU,
        _MayHaveArrayHDU,
        _MustHaveWavelengthHDU,
      ):
    _EXTNAME = 'OI_FLUX'
    
    _CARDS = [
        ('CALSTAT', True,  _u.is_calstat, None, 
            'Flux is (C)alibrated or (U)ncalibrated'),
        ('FOV',     False, _u.is_pos,     None, 
            'Field of view (arcsec)'),
        ('FOVTYPE', False, _u.is_fovtype, None, 
            'Type of field of view'),
    ]
    
    _COLUMNS = [
        ('FLUXDATA',          True,  'D',  (_NW,), None, 'any', None,
            'Flux'),
        ('FLUXERR',           True,  'D',  (_NW,), None, 'any', None,
            'Flux uncertainty'),
        ('STA_INDEX',         False, '1I', (),     None, None,  None,
            'Station index in matching OI_ARRAY table'),
        ('CORRINDX_FLUXDATA', False, '1J', (),     None, None,  None,
            'Index of 1st flux in matching OI_CORR matrix'),
    ]
    
    
    @classmethod
    def from_data(cls, *, insname, mjd, fluxdata, target_id, 
            calibrated=False,
            fov=0., fovtype='FWHM', 
            fits_keywords={}, **columns):
        """
Arguments
---------

version (int)
    version of the OIFITS standard if it cannot be determined from
    context (optional)

date (str)
    date at the start of the first observation (optional, deduced from
    mdj and int_time)
insname (str)
    name of the OI_WAVELENGTH table containing the instrumental setup 
arrname (str)         
    name of the OI_ARRAY table containing array properties (only if
    not calibrated)
corrname (str)
    name of the correlation matrix (only if uncalibrated)

mjd (float, NROWS)
    MJD of observation
int_time (float, NROWS)
    Integration time
target_id (int, NROWS)
    Target ID in the OI_TARGET table
sta_index (int, NROWS)
    station index of the flux (only if uncalibrated)

fluxdata (float, NROWS × NWAVE)
    flux
fluxerr (float, NROWS × NWAVE)
    uncertainty on the latter (optional, defaults to 0.)
flag (bool, NROWS × NWAVE)
    flag (optional, defaults to False)

corrindx_flux_data (int, NROWS)
    index in the correlation matrix for the first fluxdata (optional)

fits_keywords (dict)
    additional FITS keywords (optional)

Additional arguments
--------------------

Any additional keyword argument will be appended as a non-standard FITS 
column with its name prefixed with NS_ 

        """

        if calibrated:
            fits_keywords = dict(calstat='C', fov=fov, fovtype=fovtype,
                                 **fits_keywords) 
        else:
            fits_keywords = dict(calstat='U', **fits_keywords)
            
        columns = dict(fluxdata=fluxdata, target_id=target_id, **columns)

        return super().from_data(insname=insname, mjd=mjd,
                    fits_keywords=fits_keywords, **columns)
    
    def get_obs_type(self, name, shape='data', flatten=False):
   
        typ = f"{'un' if self.header['CALSTAT'] == 'U' else ''}calibrated flux" 
            
        return self._resize_data(typ, shape, flatten)

    def update_uv(self):
        pass
    
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

    def __mod__(self, other):

        h1, h2 = self.header, other.header

        return (super().__mod__(other) and
                h1['CALSTAT'] == h2['CALSTAT'] and
                h1.get('FOV', 0) == h2.get('FOV', 0) and
                h1.get('FOVTYPE', '') == h2.get('FOVTYPE', ''))
 

class FluxHDU1(
        _FluxHDU,
        _DataHDU21,
      ):
    """

First revision of the OI_FLUX table, OIFITS v. 2

    """ 
    pass

new_flux_hdu = _FluxHDU.from_data
