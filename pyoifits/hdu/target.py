from .table import _OITableHDU, _OITableHDU11, _OITableHDU22
from .. import utils as _u

import numpy as _np
from astroquery.simbad import Simbad as _Simbad

_deg = _np.deg2rad(1)
_milliarcsec = _deg / 3_600_000

class _MustHaveTargetHDU(_OITableHDU):

    def _get_target_field(self, name, shape='none', flatten=False,
            default=None):

        refhdu = self.get_targetHDU()
        if refhdu is None:
            val = default
        else:
            val = self._xmatch(name, refhdu, 'TARGET_ID')
        return self._resize_data(val, shape, flatten)

    def get_target(self, shape='none', flatten=False, default='N/A'):
        obj = self._container[0].header.get('OBJECT', default)
        if obj != 'MULTI':
            default = obj
        return self._get_target_field('TARGET', shape, flatten, default)
    def get_equinox(self, shape='none', flatten=False):
        return self._get_target_field('EQUINOX', shape, flatten)
    def get_ra(self, shape='none', flatten=False):
        return self._get_target_field('RAEP0', shape, flatten)
    def get_dec(self, shape='none', flatten=False):
        return self._get_target_field('DECEP0', shape, flatten)
    def get_parallax(self, shape='none', flatten=False):
        return self._get_target_field('PARALLAX', shape, flatten)
    def get_pmra(self, shape='none', flatten=False):
        return self._get_target_field('PMRA', shape, flatten)
    def get_pmdec(self, shape='none', flatten=False):
        return self._get_target_field('PMDEC', shape, flatten)
    def get_rv(self, shape='none', flatten=False):
        return self._get_target_field('SYSVEL', shape, flatten)
    def get_spectype(self, shape='none', flatten=False):
        return self._get_target_field('SPECTYP', shape, flatten)
    def get_category(self, shape='none', flatten=False):
        return self._get_target_field('CATEGORY', shape, flatten)
    
    def get_targetHDU(self):
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


def _is_veltyp(s):
    return s in ['LSR', 'HELIOCEN', 'BARYCENTR', 'TOPOCENT']

def _is_veldef(s):
    return s in ['OPTICAL', 'RADIO']

class _TargetHDU(_MustHaveTargetHDU):
    
    _EXTNAME = 'OI_TARGET' 
    _COLUMNS = [
        ('TARGET_ID',  True, '>i2',  (), _u.is_strictpos, None, None,
            'target ID for cross-reference'),
        ('RAEP0',      True, '>f8',  (), None,            None, "deg",
            'right ascension at epoch'), 
        ('DECEP0',     True, '>f8',  (), None,            None, "deg",
            'declination at epoch'), 
        ('EQUINOX',    True, '>f4',  (), None,            None, "yr",
            'equinox and epoch'),
        ('RA_ERR',     True, '>f8',  (), None,            0.,   "deg",
            'uncertainty on right ascension'),  
        ('DEC_ERR',    True, '>f8',  (), None,            0.,   "deg",
            'uncertainty on declination'),
        ('SYSVEL',     True, '>f8',  (), None,            None, "m/s",
            'radial velocity'), 
        ('VELTYP',     True, '<U8',  (), _is_veltyp,      None, None,
            'reference frame for radial velocity'), 
        ('VELDEF',     True, '<U8',  (), _is_veldef,      None, None,
            'definition for radial velocity (e.g. RADIO)'),
        ('PMRA',       True, '>f8',  (), None,            0.,   "deg/yr",
            'proper motion in right ascension'), 
        ('PMDEC',      True, '>f8',  (), None,            0.,   "deg/yr",
            'proper motion in declination'),
        ('PMRA_ERR',   True, '>f8',  (), None,            0.,   "deg/yr",
            'uncertainty on proper motion in r.a.'), 
        ('PMDEC_ERR',  True, '>f8',  (), None,            0.,   "deg/yr",
            'uncertainty on proper motion in dec.'),
        ('PARALLAX',   True, '>f4',  (), None,            None, "deg",
            'parallax'), 
        ('PARA_ERR',   True, '>f4',  (), None,            0.,   "deg",
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

    @classmethod
    def from_data(cls, *, version=2, target_id=None, target=None, 
        ra=None, dec=None, equinox=2000, ra_err=_np.nan, dec_err=_np.nan, 
        sysvel=_np.nan, veltyp='BARYCENTRIC', veldef='OPTICAL',
        pmra=None, pmdec=None, pmra_err=_np.nan, pmdec_err=_np.nan, 
        parallax=_np.nan,
        para_err=_np.nan, spectyp='UNKNOWN', category='SCI',
        fits_keywords={}, **columns):

        nrows = len(target)
        if target_id is None:
            target_id = list(range(1, nrows + 1))

        columns = dict(target=target, target_id=target_id, 
            raep0=ra, decep0=dec, equinox=equinox,
            ra_err=ra_err, dec_err=dec_err, sysvel=sysvel, veltyp=veltyp,
            veldef=veldef, pmra=pmra, pmdec=pmdec, pmra_err=pmra_err,
            pmdec_err=pmdec_err, parallax=parallax, para_err=para_err,
            spectyp=spectyp, category=category,
            **columns)

        return super().from_data(version=version, fits_keywords=fits_keywords,
                        **columns)
            
    @classmethod
    def from_simbad(cls, simbad_id, *, version=2, target_id=None,
        category='SCI', fits_keywords={}, **columns):

        simbad = _Simbad()
        simbad.remove_votable_fields(*simbad.get_votable_fields()[1:])
        simbad.add_votable_fields(
            'ra(d;A;ICRS;J2000;2000)', 'dec(d;D;ICRS;J2000;2000)',
            'coo_err_angle', 'coo_err_maja', 'coo_err_mina',
            'plx', 'plx_error',
            'rvz_radvel', 'rvz_wavelength',
            'pmra', 'pmdec', 'pm_err_angle', 'pm_err_maja', 'pm_err_mina',
            'sp',
        )
        tab = simbad.query_objects(simbad_id)

        # target ID must be ascii, if not, pick simbad
        main_id = tab['MAIN_ID']
        target = [s if ascii(s)[1:-1] == s else m.decode() 
                            for s, m in zip(simbad_id, main_id)]

        def tolist(x, deflt=_np.nan):
            x = [e.decode() if isinstance(e, bytes) else e for e in x.tolist()]
            x = _np.array([deflt if e in ['', None] else e for e in x])
            return x

        def ellipse_to_xy_err(a, b, theta):
            a2, b2 = a**2, b**2
            cos2, sin2 = _np.cos(theta) ** 2, _np.sin(theta) ** 2
            x = _np.sqrt(a2 * cos2 + b2 * sin2)
            y = _np.sqrt(a2 * sin2 + b2 * cos2)
            return x,y 

        # coordinates and parallaxes
        ra = tolist(tab['RA_d_A_ICRS_J2000_2000'])
        dec = tolist(tab['DEC_d_D_ICRS_J2000_2000'])
        a = _milliarcsec * tolist(tab['COO_ERR_MAJA'])
        b = _milliarcsec * tolist(tab['COO_ERR_MINA'])
        theta = _deg * tolist(tab['COO_ERR_ANGLE'])
        ra_err, dec_err  = ellipse_to_xy_err(a, b, theta)
        equinox = 2000.
 
        parallax = _milliarcsec * tolist(tab['PLX_VALUE'])
        para_err = _milliarcsec * tolist(tab['PLX_ERROR'])
        
        # velocity and proper motion
        sysvel = 1e3 * tolist(tab['RVZ_RADVEL'])
        veldef = tolist(tab['RVZ_WAVELENGTH'], 'OPTICAL')
        veltyp = 'BARYCENTRIC'
        
        pmra = _milliarcsec * tolist(tab['PMRA'])
        pmdec = _milliarcsec * tolist(tab['PMDEC'])
        a = _milliarcsec * tolist(tab['PM_ERR_MAJA'])
        b = _milliarcsec * tolist(tab['PM_ERR_MINA'])
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
                category=category) 

def _reshape_to_table(x, nrows):
    if not _np.shape(x):
        return _np.full(x, (nrows,))
    return _np.asarray(x)

class TargetHDU1(
        _TargetHDU,
        _OITableHDU11, # OIFITS1, table rev. 1
      ):
    _COLUMNS = [
        ('TARGET',  True, '<U16', (), _u.is_nonempty, None, None,
            'name of the celestial object'),
        ('SPECTYP', True, '<U16', (), None,           None, None,
            'spectral type'),
    ]

class TargetHDU2(
        _TargetHDU, 
        _OITableHDU22, # OIFITS2, table rev. 2
      ):
    _COLUMNS = [
        ('TARGET',   True,  '<U32', (), _u.is_nonempty, None,  None,
            'name of the celestial object'),
        ('SPECTYP',  True,  '<U32', (), None,           None,  None,
            'spectral type'),
        ('CATEGORY', False, '<U3',  (), _u.is_category, 'SCI', None,
            'observation category: SCIence or CALibration'),
    ]

new_target_hdu = _TargetHDU.from_data
new_target_hdu_from_simbad = _TargetHDU.from_simbad
