"""
Implementation of the OI_TARGET binary table extension
"""

from .table import _OITableHDU, _OITableHDU11, _OITableHDU22
from .. import utils as _u

import numpy as _np
from astroquery.simbad import Simbad as _Simbad
from astropy.coordinates import SkyCoord as _SkyCoord
from astropy import units as _units
from astropy.time import Time as _Time

_deg = _np.deg2rad(1)
_arcsec = _deg / 3600
_milliarcsec = _arcsec / 1000

def _decode(s):
    if isinstance(s, bytes):
        return s.decode()
    return s

class _MustHaveTargetHDU(_OITableHDU):

    def _get_target_field(self, name, shape='none', flatten=False,
            default=None):

        refhdu = self.get_targetHDU()
        if refhdu is None:
            val = default
        else:
            val = self._xmatch(refhdu, 'TARGET_ID', name=name)
        return self._resize_data(val, shape, flatten)

    def get_target(self, shape='none', flatten=False, default='N/A'):
        """

Get the target name associated to each row of the table

Arguments
---------

shape (str)
    Shape of the returned array
    * None or 'table':  1D array (str × NROWS), one per row.
    * 'data': 2D array (str × NOBS × NWAVE), matching the 
        interferometric data for OI_VIS, OI_VIS2, OI_T3, and 
        OI_FLUX tables; 1D array otherwise.

flatten (bool)
    Flattens to 1D array.

        """
        obj = self._container[0].header.get('OBJECT', default)
        if obj != 'MULTI':
            default = obj
        return self._get_target_field('TARGET', shape, flatten, default)

    def get_sky_coord(self):
        """
Get the full coordinates of the targets, including distance and 
motions.

Returns
-------

astropy.coordinates.SkyCoord object

        """
        coordinates = []

        refhdu = self.get_targetHDU()
        indices = self._xmatch(refhdu, 'TARGET_ID')
        rows = refhdu.data[indices]        

        # A bug in the OIFITS standard (explicit in v. 1)doesn't make a
        # distinction between epoch and equinox and impose both to be a given
        # value.  This excludes ICRS where equinox is fixed.  With the
        # recommended value 2000.0, it means FK5 J2000.0/2000.0 ≃ ICRS with ~70
        # mas accuracy.  Not interferometry-grade precision...
        equinox = _Time(rows['EQUINOX'], format='decimalyear')
        epoch = equinox
        frame = 'fk5'

        ra = rows['RAEP0'] * _units.deg
        dec = rows['DECEP0'] * _units.deg
        pm_ra = rows['PMRA'] * _units.deg / _units.yr
        pm_dec = rows['PMDEC'] * _units.deg / _units.yr
        rv = rows['SYSVEL'] * _units.m / _units.s
        plx = rows['PARALLAX'] * _units.deg / _units.yr
        plx[plx == 0] = _np.nan

        dist = 1 / plx.to_value('arcsec / yr') * _units.pc
        
        coo = _SkyCoord(
            ra=ra, dec=dec, distance=dist,
            pm_ra_cosdec=pm_ra, pm_dec=pm_dec, radial_velocity=rv,
            equinox=equinox, obstime=epoch, frame=frame,
        )


        return coo
        
    def get_equinox(self, shape='none', flatten=False):
        """

Get the equinox and epoch associated to each row of the table.

Arguments
---------

shape (str)
    Shape of the returned array
    * None or 'table':  1D array (str × NROWS), one per row.
    * 'data': 2D array (str × NOBS × NWAVE), matching the 
        interferometric data for OI_VIS, OI_VIS2, OI_T3, and 
        OI_FLUX tables; 1D array otherwise.

flatten (bool)
    Flattens to 1D array.

        """
        return self._get_target_field('EQUINOX', shape, flatten)
    def get_ra(self, shape='none', flatten=False):
        """

Get the right ascension associated to each row of the table

Arguments
---------

shape (str)
    Shape of the returned array
    * None or 'table':  1D array (str × NROWS), one per row.
    * 'data': 2D array (str × NOBS × NWAVE), matching the 
        interferometric data for OI_VIS, OI_VIS2, OI_T3, and 
        OI_FLUX tables; 1D array otherwise.

flatten (bool)
    Flattens to 1D array.

        """
        return self._get_target_field('RAEP0', shape, flatten)
    def get_dec(self, shape='none', flatten=False):
        """

Get the declination associated to each row of the table

Arguments
---------

shape (str)
    Shape of the returned array
    * None or 'table':  1D array (str × NROWS), one per row.
    * 'data': 2D array (str × NOBS × NWAVE), matching the 
        interferometric data for OI_VIS, OI_VIS2, OI_T3, and 
        OI_FLUX tables; 1D array otherwise.

flatten (bool)
    Flattens to 1D array.

        """
        return self._get_target_field('DECEP0', shape, flatten)
    def get_parallax(self, shape='none', flatten=False):
        """

Get the parallax associated to each row of the table

Arguments
---------

shape (str)
    Shape of the returned array
    * None or 'table':  1D array (str × NROWS), one per row.
    * 'data': 2D array (str × NOBS × NWAVE), matching the 
        interferometric data for OI_VIS, OI_VIS2, OI_T3, and 
        OI_FLUX tables; 1D array otherwise.

flatten (bool)
    Flattens to 1D array.

        """
        return self._get_target_field('PARALLAX', shape, flatten)
    def get_pmra(self, shape='none', flatten=False):
        """

Get the proper motion in right ascension associated to each row of 
the table

Arguments
---------

shape (str)
    Shape of the returned array
    * None or 'table':  1D array (str × NROWS), one per row.
    * 'data': 2D array (str × NOBS × NWAVE), matching the 
        interferometric data for OI_VIS, OI_VIS2, OI_T3, and 
        OI_FLUX tables; 1D array otherwise.

flatten (bool)
    Flattens to 1D array.

        """
        return self._get_target_field('PMRA', shape, flatten)
    def get_pmdec(self, shape='none', flatten=False):
        """

Get the proper motion in declination associated to each row of the table

Arguments
---------

shape (str)
    Shape of the returned array
    * None or 'table':  1D array (str × NROWS), one per row.
    * 'data': 2D array (str × NOBS × NWAVE), matching the 
        interferometric data for OI_VIS, OI_VIS2, OI_T3, and 
        OI_FLUX tables; 1D array otherwise.

flatten (bool)
    Flattens to 1D array.

        """
        return self._get_target_field('PMDEC', shape, flatten)
    def get_rv(self, shape='none', flatten=False):
        """

Get the radial velocity associated to each row of the table

Arguments
---------

shape (str)
    Shape of the returned array
    * None or 'table':  1D array (str × NROWS), one per row.
    * 'data': 2D array (str × NOBS × NWAVE), matching the 
        interferometric data for OI_VIS, OI_VIS2, OI_T3, and 
        OI_FLUX tables; 1D array otherwise.

flatten (bool)
    Flattens to 1D array.

        """
        return self._get_target_field('SYSVEL', shape, flatten)
    def get_spectype(self, shape='none', flatten=False):
        """

Get the radial velocity reference (helocientric, barycentric, etc.)
associated to each row of the table

Arguments
---------

shape (str)
    Shape of the returned array
    * None or 'table':  1D array (str × NROWS), one per row.
    * 'data': 2D array (str × NOBS × NWAVE), matching the 
        interferometric data for OI_VIS, OI_VIS2, OI_T3, and 
        OI_FLUX tables; 1D array otherwise.

flatten (bool)
    Flattens to 1D array.

        """
        return self._get_target_field('SPECTYP', shape, flatten)
    def get_category(self, shape='none', flatten=False):
        """

Get the target category (SCI/CAL) associated to each row of the table

Arguments
---------

shape (str)
    Shape of the returned array
    * None or 'table':  1D array (str × NROWS), one per row.
    * 'data': 2D array (str × NOBS × NWAVE), matching the 
        interferometric data for OI_VIS, OI_VIS2, OI_T3, and 
        OI_FLUX tables; 1D array otherwise.

flatten (bool)
    Flattens to 1D array.

        """
        return self._get_target_field('CATEGORY', shape, flatten)
    
    def get_targetHDU(self):
        """

Get the corresponding OI_TARGET HDU. 

        """
        return self._container.get_targetHDU()

    def _verify(self, option='warn'):

        errors = super()._verify(option)

        # Verify Target ID is correct (> 0) and referenced

        val = _np.unique(self.TARGET_ID)
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

class _TargetHDU(_MustHaveTargetHDU):
    
    _EXTNAME = 'OI_TARGET' 
    _COLUMNS = [
        ('TARGET_ID',  True, 'I',  (), _u.is_strictpos, None, None,
            'target ID for cross-reference'),
        ('RAEP0',      True, '1D',  (), None,            None, "deg",
            'right ascension at epoch'), 
        ('DECEP0',     True, '1D',  (), None,            None, "deg",
            'declination at epoch'), 
        ('EQUINOX',    True, '1E',  (), None,            None, "yr",
            'equinox and epoch'),
        ('RA_ERR',     True, '1D',  (), None,            0.,   "deg",
            'uncertainty on right ascension'),  
        ('DEC_ERR',    True, '1D',  (), None,            0.,   "deg",
            'uncertainty on declination'),
        ('SYSVEL',     True, '1D',  (), None,            None, "m/s",
            'radial velocity'), 
        ('VELTYP',     True, '8A',  (), _u.is_veltyp,    None, None,
            'reference frame for radial velocity'), 
        ('VELDEF',     True, '8A',  (), _u.is_veldef,    None, None,
            'definition for radial velocity (e.g. RADIO)'),
        ('PMRA',       True, '1D',  (), None,            0.,   "deg/yr",
            'proper motion in right ascension'), 
        ('PMDEC',      True, '1D',  (), None,            0.,   "deg/yr",
            'proper motion in declination'),
        ('PMRA_ERR',   True, '1D',  (), None,            0.,   "deg/yr",
            'uncertainty on proper motion in r.a.'), 
        ('PMDEC_ERR',  True, '1D',  (), None,            0.,   "deg/yr",
            'uncertainty on proper motion in dec.'),
        ('PARALLAX',   True, '1E',  (), None,            None, "deg",
            'parallax'), 
        ('PARA_ERR',   True, '1E',  (), None,            0.,   "deg",
            'error on parallax'),
    ]
    
    ang_dist_max = 5e-9 # approx 1 mas
    
    def _verify(self, option='warn'):

        errors = super()._verify(option)

        target_id = self.TARGET_ID
        if len(_np.unique(target_id)) == len(target_id):
            return errors

        err_text = f"Repeated TARGET_ID in {type(self).__name__}"
        err = self.run_option(option, err_text, fixable=False)
        errors.append(err)

        return errors
   
    def __add__(self, other):

        return self._merge_helper(other)

    def __mod__(self, other):
        
        return self & other
               
    def is_referred_to_by(self, other):
        return (not isinstance(other, _TargetHDU) and
                isinstance(other, _MustHaveTargetHDU))

    def merge(self, *others):
   
        container = self.get_container()
        dist_max = container._merge_target_distance
        name_match = container._merge_target_name_match
 
        def eq(x, y):
            return ((not name_match or x['TARGET'] == y['TARGET']) and
                    (dist_max > 360 or x['EQUINOX'] == y['EQUINOX']) and
                    abs(x['RAEP0'] - y['RAEP0']) <= dist_max and
                    abs(x['DECEP0'] - y['DECEP0']) <= dist_max)
                 
        return self._merge_helper(*others, id_name='TARGET_ID', equality=eq)

    def _trim_helper(self, *, target_filter=lambda targ: True, 
            wave_filter=None, insname_filter=None, keep_ns_columns=False):

        tkeep = _np.vectorize(target_filter)(self.get_target())
        
        columns = {}
        data = self.data
        standard_colnames = self._get_oi_colnames()
        for column in self.columns:
            colname = column.name
            if not keep_ns_columns and colname not in standard_colnames:
                continue
            columns[colname.lower()] = data[colname][tkeep]

        thdu = self._from_data(fits_keywords=self.header, **columns)

        return thdu



    @classmethod
    def from_data(cls, *, target_id=None, target, ra, dec,  
            fits_keywords={}, **columns):
        """

Build an OI_TARGET extension from data.  In the following NTARGET is
the number of targets.

Arguments
---------

version (int, optional, default: 2)
    Version of the OIFITS standard to use if it cannot be determined
    from context

target (str × NTARGET):
    Target name.
target_id (int × NTARGET, optional, defaults to 1 .. NTARGET)
    Target ID for cross-reference. 
ra (float × NTARGET)
    Right ascension in degrees
ra_err (float × NTARGET, optional, default: NaN)
    Uncertainty on the latter
dec (float × NTARGET)
    Declination in degrees
dec_err (float × NTARGET, optional, default: NaN)
    Uncertainty on the latter
parallax (float × NTARGET, optional, default: NaN)
    Parallax in degrees
para_err (float × NTARGET, optional, default: NaN)
    Uncertainty on the latter
equinox (float × NTARGET, optional, default: 2000)
    Equinox and epoch of observation in the FK5 frame.
pmra (float × NTARGET, optional, default: NaN)
    Proper movement in right ascension in degrees per yer.
pmra_err (float × NTARGET, optional, default: NaN)
    Uncertainty on the latter
pmdec (float × NTARGET, optional, default: NaN)
    Proper movement in declination in degrees
dec_err (float × NTARGET, optional, default: NaN)
    Uncertainty on the latter
sysvel (float × NTARGET, optional, default: NaN)
    Radial velocity
veltyp (str  × NTARGET, optional, default: 'BARYCENTRIC')
    Radial velocity reference frame
veldef (str  × NTARGET, optional, default: 'OPTICAL')
    Wavelength at which velocity was obtained
spectyp (str  × NTARGET, optional, default: 'UNKNOWN')
    Spectral type

Arguments in OIFITS2 only
-------------------------

category (str  × NTARGET, optional)
    Target category (SCI: science, CAL: calibrator)

Additional arguments
--------------------

Any additional keyword argument will be appended as a non-standard FITS 
column with its name prefixed with NS_ 

Warning
-------

The standard mandates EQUINOX be equal to both equinox and epoch
(see version 1, version does not modify any of it). It rules out the 
ICRS.

1. For applications needing interferometric-grade precision of the 
coordinates, take care that the FK5 frame is used, not ICRS, with
~70 milliarcseconds differences.

2. FK5 J2000.0 or ICRS coordinates with a different epoch (e.g. GAIA)
cannot be entered directly and must be converted to an FK5 system
where equinox and epoch are equal, preferably 2000.0

        """
        nrows = len(target)
        if target_id is None:
            target_id = list(range(1, nrows + 1))

        _u.store_default(columns, 'ra_err', default=_np.nan)
        _u.store_default(columns, 'dec_err', default=_np.nan)
        _u.store_default(columns, 'pmra', default=_np.nan)
        _u.store_default(columns, 'pmdec', default=_np.nan)
        _u.store_default(columns, 'pmra_err', default=_np.nan)
        _u.store_default(columns, 'pmdec_err', default=_np.nan)
        _u.store_default(columns, 'parallax', default=_np.nan)
        _u.store_default(columns, 'para_err', default=_np.nan)
        _u.store_default(columns, 'sysvel', default=_np.nan)
        _u.store_default(columns, 'veltyp', default='BARYCENTRIC')
        _u.store_default(columns, 'veldef', default='OPTICAL')
        _u.store_default(columns, 'equinox', default=2000.0)

        columns = dict(target=target, raep0=ra, decep0=dec,
                target_id=target_id, **columns) 

        return super().from_data(fits_keywords=fits_keywords, **columns)
            
    @classmethod
    def from_simbad(cls, simbad_id, *, version=2, target_id=None,
        fits_keywords={}, **columns):
        """

Create an OI_TARGET table from simbad identifiers

Arguments
---------

version (int, optional)
    Version of the OIFITS standard if it cannot be deduced from context
simbad_id (str × NTARGET)
    IDs recognised by SIMBAD
target_id (int × NTARGET, optional, default: 1 .. NTARGET)
    target ID for cross-reference 
category (str × NTARGET, optional)
    target catogory, SCI or CAL
fits_keywords (dict, optional, default: {})
    additional FITS header keywords 

Additional arguments
--------------------

Any additional keyword argument will be appended as a non-standard FITS 
column with its name prefixed with NS_ 

        """

        simbad = _Simbad()
        simbad.remove_votable_fields(*simbad.get_votable_fields()[1:])
        # OIFITS standard imposes equinox = epoch for the frame which
        # de facto excludes ICRS, we use FK5.
        simbad.add_votable_fields(
            'ra(d;A;FK5;J2000;2000)', 'dec(d;D;FK5;J2000;2000)',
            'coo_err_angle', 'coo_err_maja', 'coo_err_mina',
            'plx', 'plx_error',
            'rvz_radvel', 'rvz_wavelength',
            'pmra', 'pmdec', 'pm_err_angle', 'pm_err_maja', 'pm_err_mina',
            'sp',
        )
        tab = simbad.query_objects(simbad_id)

        # target ID must be ascii, if not, pick simbad
        main_id = tab['MAIN_ID']
        target = [s if ascii(s)[1:-1] == s else _decode(m)
                            for s, m in zip(simbad_id, main_id)]

        def tolist(x, deflt=_np.nan):
            x = [_decode(e) for e in x.tolist()]
            x = _np.array([deflt if e in ['', None] else e for e in x])
            return x

        def ellipse_to_xy_err(a, b, theta):
            a2, b2 = a**2, b**2
            cos2, sin2 = _np.cos(theta) ** 2, _np.sin(theta) ** 2
            x = _np.sqrt(a2 * cos2 + b2 * sin2)
            y = _np.sqrt(a2 * sin2 + b2 * cos2)
            return x,y 

        # coordinates and parallaxes
        ra = tolist(tab['RA_d_A_FK5_J2000_2000'])
        dec = tolist(tab['DEC_d_D_FK5_J2000_2000'])
        a = _milliarcsec / _deg * tolist(tab['COO_ERR_MAJA'])
        b = _milliarcsec  / _deg * tolist(tab['COO_ERR_MINA'])
        theta = _deg * tolist(tab['COO_ERR_ANGLE'])
        ra_err, dec_err  = ellipse_to_xy_err(a, b, theta)
        equinox = 2000.
 
        parallax = _milliarcsec / _deg * tolist(tab['PLX_VALUE'])
        para_err = _milliarcsec / _deg * tolist(tab['PLX_ERROR'])
        
        # velocity and proper motion
        sysvel = 1e3 * tolist(tab['RVZ_RADVEL'])
        veldef = tolist(tab['RVZ_WAVELENGTH'], 'OPTICAL')
        veltyp = 'BARYCENTRIC'
        
        pmra = _milliarcsec / _deg * tolist(tab['PMRA'])
        pmdec = _milliarcsec / _deg * tolist(tab['PMDEC'])
        a = _milliarcsec / _deg * tolist(tab['PM_ERR_MAJA'])
        b = _milliarcsec / _deg * tolist(tab['PM_ERR_MINA'])
        theta = _deg * tolist(tab['PM_ERR_ANGLE'])
        pmra_err, pmdec_err = ellipse_to_xy_err(a, b, theta)
       
        # spectral type
        spectyp = tolist(tab['SP_TYPE'], 'UNKNOWN') 

        return cls.from_data(version=version, target_id=target_id,
                target=target, ra=ra, dec=dec, equinox=equinox,
                ra_err=ra_err, dec_err=dec_err,
                sysvel=sysvel, veltyp=veltyp, veldef=veldef,
                pmra=pmra, pmdec=pmdec, pmra_err=pmra_err, pmdec_err=pmdec_err,
                parallax=parallax, para_err=para_err, spectyp=spectyp,
                **columns)

def _reshape_to_table(x, nrows):
    if not _np.shape(x):
        return _np.full(x, (nrows,))
    return _np.asarray(x)

class TargetHDU1(
        _TargetHDU,
        _OITableHDU11, # OIFITS1, table rev. 1
      ):
    """

First revision of the OI_TARGET binary table, OIFITS v. 1

    """
    _COLUMNS = [
        ('TARGET',  True, '16A', (), _u.is_nonempty, None, None,
            'name of the celestial object'),
        ('SPECTYP', True, '16A', (), None,           None, None,
            'spectral type'),
    ]

class TargetHDU2(
        _TargetHDU, 
        _OITableHDU22, # OIFITS2, table rev. 2
      ):
    """

Second revision of the OI_TARGET binary table, OIFITS v. 2

    """
    _COLUMNS = [
        ('TARGET',   True,  '32A', (), _u.is_nonempty, None,  None,
            'name of the celestial object'),
        ('SPECTYP',  True,  '32A', (), None,           None,  None,
            'spectral type'),
        ('CATEGORY', False, '3A',  (), _u.is_category, 'SCI', None,
            'observation category: SCIence or CALibration'),
    ]

new_target_hdu = _TargetHDU.from_data
new_target_hdu_from_simbad = _TargetHDU.from_simbad
