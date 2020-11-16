from .data import _DataHDU11, _DataHDU22
from .t2 import _T2HDU
from .wavelength import _NW
from .. import utils as _u

import numpy as np

class _VisHDU(_T2HDU):
    """
OI_VIS binary table extension containing the (possibly differential) 
visibility amplitudes, phases, as well as correlated fluxes.
    """
    _EXTNAME = 'OI_VIS'    
    _COLUMNS = [
        ('VISAMP',    True, 'D', (_NW,), None, None, None,
            'visibility amplitude'), 
        ('VISAMPERR', True, 'D', (_NW,), None, None, None,
            'uncertainty on visibility amplitude'),
        ('VISPHI',    True, 'D', (_NW,), None, None, 'deg',
            'phase'), 
        ('VISPHIERR', True, 'D', (_NW,), None, None, 'deg', 
            'uncertainty on phase'),
    ]
    
    @classmethod
    def from_data(cls, *, visamp, visphi, fits_keywords={}, **columns):
        """

Build an OI_VIS table from data.  In the following NWAVE indicates the
number of spectral channels in the configuration and NOBS the number
of observations.

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
    name of the OI_ARRAY table containing array properties (optional in
    version 1)
fits_keywords (dict)
    additional FITS keywords (optional)

mjd (float, NOBS)
    MJD of observation
int_time (float, NOBS)
    Integration time
target_id (int, NOBS)
    Target ID in the OI_TARGET table
sta_index (int, NOBS × 2)
    Station IDs of the telescope pair in the OI_ARRAY target (optional in
    version 1)
ucoord, vcoord(float, NOBS)
    (u, v) coordinates of the  baselines

visamp (float, NOBS × NWAVE)
    visibility amplitudes
visamperr (float, NOBS × NWAVE)
    uncertainty on the latter (optional, defaults to 0.)
visphi (float, NOBS × NWAVE)
    visibility phases
visphi (float, NOBS × NWAVE)
    visibility phases (optional, defaults to 0.)
flag (bool, NOBS × NWAVE)
    flag for data (optional, defaults to False)

Arguments in OIFITS1 only
-------------------------

time (float, NOBS)
    seconds since midnight UT on first observation (DATE-OBS i.e. date)

Arguments in OIFITS2 only
-------------------------

corrname (str)
    name of the correlation matrix (optional)
amptyp (str)
    type of amplitude: absolute, differential, correlated flux (optional)
phityp (str)
    type of phase: absolute, differential (optional)
amporder (int)
    order of the polynomial fit for differential amplitudes (optional)
phiorder (int)
    order of the polynomial fit for differential phases (optional)

rvis (float, NOBS × NWAVE)
    real part of the coherent flux (optional)
rviserr (float, NOBS × NWAVE)
    uncertainty on the latter (optional, defaults to 0.)
ivis (float, NOBS × NWAVE)
    imaginary part of the coherent flux (optional)
iviserr (float, NOBS × NWAVE)
    uncertainty on the latter (optional, defaults to 0.)
corrindx_visamp (int, NOBS × NWAVE)
    index in the correlation matrix for the first visamp (optional)  
corrindx_visphi (int, NOBS × NWAVE)
    index in the correlation matrix for the first visphi (optional)
corrindx_rvis (int, NOBS × NWAVE)
    index in the correlation matrix for the first rvis (optional)
corrindx_ivis (int, NOBS × NWAVE)
    index in the correlation matrix for the first ivis (optional)

visref_map (int, NOBS × NWAVE × NWAVE)
    indices of the spectral channels used to determine the differential
    quantities

Additional arguments
--------------------

Any additional keyword argument will be appended as a non-standard FITS 
column with its name prefixed with NS_ 

        """
        columns = dict(visamp=visamp, visphi=visphi,  **columns)

        shape = self._get_columns_shape(**columns)
        for name in ['visamp', 'visphi', 'rvis', 'ivis']:
            if name in columns:
                _u.store_default(columns, f"{name}err", 0., shape)

        return super().from_data(fits_keywords=fits_keywords, **columns)

class VisHDU1(
        _VisHDU,
        _DataHDU11,
      ):
    """

First revision of the OI_VIS binary table, OIFITS v. 1.

    """
    pass

_spos = _u.is_strictpos

class VisHDU2(
        _VisHDU,
        _DataHDU22,
      ):
    """

Second revision of the OI_VIS binary table, OIFITS v. 2.

    """
    _COLUMNS = [
        ('CORRINDX_VISAMP',  False, '1J', (),        _spos, None, None,
            'index on 1st amp. in matching OI_CORR matrix'), 
        ('CORRINDEX_VISPHI', False, '1J', (),        _spos, None, None,
            'index on 1st phase in matching OI_CORR matrix'),
        ('RVIS',             False, 'D', (_NW,),     None,  None, 'any', 
            'real part of correlated flux'),                     
        ('RVISERR',          False, 'D', (_NW,),     None,  None, 'any', 
            'uncertainty on real part of correlated flux'),      
        ('IVIS',             False, 'D', (_NW,),     None,  None, 'any', 
            'imag. part of correlated flux'),                    
        ('IVISERR',          False, 'D', (_NW,),     None,  None, 'any', 
            'uncertainty on imag. part of correlated flux'),
        ('CORRINDX_RVIS',    False, '1J', (),        _spos, None, None,
            'index of 1st RVIS in matching OI_CORR matrix'),
        ('CORRINDEX_IVIS',   False, '1J', (),        _spos, None, None,
            'index of 1st IVIS in matching OI_CORR matrix'),
        ('VISREFMAP',        False, 'L', (_NW,_NW,), None,  None, None,
            'reference channels for differential quantities'),
    ]
    
    _CARDS = [
        ('AMPTYP',   False, _u.is_amptyp2,  None, 
            'amplitude type'),
        ('PHITYP',   False, _u.is_phityp2,  None, 
            'phase is absolute or differential'),
        ('AMPORDER', False, _u.is_posint, 0, 
            'polynomial order in diff. amplitude'),
        ('PHIORDER', False, _u.is_posint, 0, 
            'polynomial order in diff. phase'),
    ]
    
    def get_reference_channel(self, shape='data', flatten=False):

        visrefmap = getattr(self, 'VISREFMAP', 0)
        if visrefmap == 0:
            return super().get_reference_channel('data', False)

        nwave = visrefmap.shape[-1]
        visrefmap = np.array(visrefmap, dtype=object) # allows long int's
        visrefmap = sum(visrefmap[:,:,i] << (i + 1) for i in range(nwave))
        if flatten:
            visrefmap = visrefmap.ravel()
 
        return visrefmap
    
    def get_obs_type(self, name, shape='data', flatten=False):
        
        if name == 'VISAMP':
            typ = self.header.get('AMPTYP', 'N/A')
        elif name == 'VISPHI':
            typ = self.header.get('PHITYP', 'N/A')
        else:
            typ = 'correlated flux'
        return self._resize_data(typ, shape, flatten)
    
    @classmethod
    def from_data(cls, *, visamp, visphi,  

        amptyp=None, phityp=None, amporder=None, phiorder=None, 
        fits_keywords={}, **columns):

        fits_keywords = dict(amptyp=amptyp, phityp=phityp, amporder=amporder,
            phiorder=phiorder, **fits_keywords)
        return super().from_data(cls, visamp=visamp, visphi=visphi, 
                        fits_keywords=fits_keywords, **columns)

    def _verify(self, option='warn'):

        errors = super()._verify(option=option)

        # Format doesn't have consistent units for VISAMP & VISAMPERR
        # Deal with it manually

        amptyp = self.header.get('AMPTYP', 'absolute')
        
        for name in ['VISAMP', 'VISAMPERR']:
            column = self.columns[name]
            if amptyp == 'correlated flux':
                if column.unit is None:
                    err_txt = f"{name} must have unit for {amptyp}"
                    err = self.run_option(option, err_txt)
                    errors.append(err)
            else:
                if column.unit not in ['', None]:
                    err_txt = f"{name} cannot have unit for AMPTYP='{amptyp}'"
                    fix_txt = "Unit removed"
                    def fix(col=column): col.unit = None
                    err = self.run_option(option, err_txt, fix_txt, fix)
                    errors.append(err)
        
        return errors

    def __mod__(self, other):

        h1, h2 = self.header, other.header

        return (super().__mod__(other) and
                h1.get('AMPTYP', '') == h2.get('AMPTYP', '') and
                h1.get('PHITYP', '') == h2.get('PHITYP', '') and
                h1.get('AMPORDER', 0) == h2.get('AMPORDER', 0) and
                h1.get('PHIORDER', 0) == h2.get('PHIORDER', 0))
        
new_vis_hdu = _VisHDU.from_data
