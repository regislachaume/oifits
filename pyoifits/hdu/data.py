from .table import _OITableHDU, _OITableHDU1, _OITableHDU2
from .table import _OIFITS1HDU, _OIFITS2HDU
from .array import _MayHaveArrayHDU,_MustHaveArrayHDU
from .target import _MustHaveTargetHDU
from .wavelength import _MustHaveWavelengthHDU
from .corr import _MayHaveCorrHDU
from .inspol import _MayHaveInspolHDU

from .. import utils as _u
from .. import coo as _coo

from astropy import table as _table
from astropy import units as _units
from numpy import ma as _ma
import re as _re
import numpy as _np
from astropy.time import Time as _Time
from scipy.spatial.transform import Rotation as _rotation
from astropy.coordinates import SkyCoord as _SkyCoord

class _DataHDU(_OITableHDU):
    
    _CARDS = [
        ('DATE-OBS', True, _u.is_nonempty, None, 'Date at start of observation'),
    ]
    _COLUMNS = [
        ('TARGET_ID', True, '1I', (),         _u.is_strictpos, None,  None,
            'Target ID for cross-reference'), 
        ('MJD',       True, '1D', (),         None,            None,  "d",
            'Modified Julian Day at start of observation'),
        ('INT_TIME',  True, '1D', (),         None,            None,  "s",
            'Integration time'), 
        ('FLAG',      True, 'L',  ('NWAVE',), None,            False, None,
            'Flag for bad quality'),
    ]
    
    def get_obs_type(self, name, shape='data', flatten=False):
        """

Get the type of observable.

Arguments
---------

name (str):
    Name of the observable (VISAMP, VISPHI, IVIS, RVIS, FLUXDATA, VIS2DATA,
    T3AMP, T3PHI).

shape (str):
    Shape of the returned argument (optional):
        * 'none': returns a scalar
        * 'table': returns a 1D array with the same length as the table 
            (NROWS)
        * 'data': returns an array with the same shape as the observable 
            (NROWS × NWAVE) 

Returns
-------

Type of observable (uncalibrated flux, calibrated flux, correlated flux, 
absolute, differential).

        """
        return self._resize_data('N/A', shape, flatten)

    def get_uv(self, shape='table', flatten=False):
        """

Get the (u, v) coordinates.

Arguments
---------

shape (str):
    Shape of the returned argument (optional):
        * 'table': returns a 2D array with th
            (2 × NROWS or 4 × NROWS)
        * 'data': returns an array with the same shape as the observable 
            (2 × NROWS × NWAVE or 4 × NROWS × NWAVE) 

Returns
-------

uv (array of float):
    (u, v) for 2T observables and (u1, v1, u2, v2) for 3T observables.


        """
        uv = [self._resize_data(self.data[x], shape_flatten)
                    for x in self._COORD_COLUMNS]

        return _np.array(uv)

    def _update_targetid(self, index_map):

        hdu = _np.copy(self)
        hdu.TARGET_ID = [index_map[id] for id in hdu.TARGET_ID]
        return hdu
    
    def _table_colnames(self, full_uv=False): 

        COORD = self._get_uvcoord_names(full_uv=full_uv)
            
        return [
            'TARGET', 'CHANNEL', 'REF_CHANNEL_BITFIELD',
            'EFF_WAVE', 'EFF_BAND', *COORD, 
            'observable', 'type', 'value', 'error', 'INSNAME',
            'ARRNAME', 'STA_CONFIG', 'MJD', 'INT_TIME'
        ]

    def _table_cols(self, full_uv=False):
        
        names = self._table_colnames(full_uv=full_uv)
        cols = []

        colnames = self.columns.names
        obs_names = [n for n in self.get_observable_names() if n in colnames]
        err_names = [n for n in self.get_error_names() if n in colnames]

        def getf(n): return self.get_field(n, 'data', True, default=0)
        def gett(n): return self.get_obs_type(n, 'data', True)
        def resize(x): return self._resize_data(x, 'data', True)
        def hstack(x): return _ma.hstack(x)
        
        for name in names:
            if name == 'value':
                col = hstack([getf(n) for n in obs_names])
            elif name == 'error':
                col = hstack([getf(n) for n in err_names])
            elif name == 'observable':
                col = _np.hstack([resize(n) for n in obs_names])
            elif name == 'type':
                col = _np.hstack([gett(n) for n in obs_names])
            else:
                col = hstack([getf(name)] * len(obs_names))
            cols.append(col)

        return cols

    def get_field(self, name, shape='none', flatten=False, default=None):

        if name == 'INSNAME':
            return self.get_insname(shape, flatten)
        if name == 'ARRNAME':
            return self.get_arrname(shape, flatten, default)
        if name == 'CORRNAME':
            return self.get_corrname(shape, flatten, default)
        if name == 'TARGET':
            return self.get_target(shape, flatten, default)
        if name == 'EFF_WAVE':
            return self.get_wave(shape, flatten)
        if name == 'EFF_BAND':
            return self.get_band(shape, flatten)
        if name == 'CHANNEL':
            return self.get_channel(shape, flatten)
        if name == 'REF_CHANNEL_BITFIELD':
            return self.get_reference_channels(shape, flatten)
        if name == 'STA_CONFIG':
            return self.get_sta_config(shape, flatten)

        DATACOLS = self._get_spec_colnames()
        DATACOLS.remove('FLAG')

        if name not in DATACOLS:
            if hasattr(self, name):
                x = getattr(self, name)
            else:
                x = default
            x = self._resize_data(x, shape, flatten)
            x = _ma.masked_array(x, mask=not hasattr(self, name))
            return x

        mask = self.FLAG
        if hasattr(self, name):
            x = getattr(self, name)
        else:
            x = self._resize_data(default, shape)
            flag = flag | True
        x = _ma.masked_array(x, mask=mask)
        if flatten:
            x = x.ravel()
 
        return x 

    def get_reference_channels(self, shape='data', flatten=False):
        """

Get the reference channels for differential quantities.  

Arguments:
----------

shape (str):
    Shape of the returned data (optional, defaults to 'data').
    * 'data': shape is NOBS × NWAVE 
    * 'table': shape NOBS

flatten (bool, optional, default: False):
   
    Whether result should be flatten to 1D array.

Returns:
--------

A bitfield (int) indicated which wavelength channel numbers are included.  
If the quantity is not differential, 0 is returned.

        """

        visref = self._resize_data(0, 'data', flatten)
        return _ma.masked_array(visref, mask=True)

    def _to_table(self, full_uv=False):

        names = self._table_colnames(full_uv=full_uv)
        cols = self._table_cols(full_uv=full_uv)
        return _table.Table(cols, names=names)

    @classmethod
    def _get_uvcoord_names(cls, full_uv=False):
       
        if full_uv:
            return ['U1COORD', 'V1COORD', 'U2COORD', 'V2COORD']
         
        names = [c for c in cls._COLUMNS['name'] if _re.match('[UV].?COORD', c)]
        names = sorted(names, key=lambda x: x[1::-1])
        return names

    def to_table(self):
        """

Convert to a flat table.

Arguments:
----------
    
full_uv (bool):
    Whether to systematically include (u1, v1, u2, v2) coordinates even
    if some of them are not relevant (e.g. flux or 2T data).  In that case
    they are NaN.

Returns:
--------

    An astropy.table.Table with one scalar observable per line.

        """
        tab = self._to_table()
        coord_names = self._get_uvcoord_names()

        for x in ['INT_TIME', *coord_names]:
            tab.columns[x].format = '7.3f'
        for x in ['EFF_WAVE', 'EFF_BAND']:
            tab.columns[x].format = '7.5e'
        tab.columns['MJD'].format = '7.5f'
        for x in ['value', 'error']:
            tab.columns[x].format = '7.5g'
        
        return tab

    def merge(self, *others):
        """

        Merge a set of OI extensions

        Example: hdu = hdu1.merge(hdu2, hdu3)

        """
        for ref in ['ARRNAME', 'INSNAME', 'CORRNAME', 'EXTNAME']:
            sref = self.header.get(ref, '')
            for other in others:
                oref = other.header.get(ref, '')
                if sref != oref:
                    txt = f'Cannot merge tables with {ref} = {sref} & {oref}'
                    raise RuntimeError(txt)

        return self._merge_helper(*others)

    @classmethod
    def from_data(cls, *, insname, arrname=None, corrname=None,
        date=None, mjd, int_time=0., flag=False, fits_keywords={}, **columns):

        fits_keywords = dict(arrname=arrname, insname=insname, 
                             corrname=corrname, **fits_keywords)
        columns = dict(mjd=mjd, int_time=int_time, flag=flag, **columns)

        if date is None:
            mjdobs = min(mjd - int_time / 86400. / 2)
            date = _Time(mjdobs, format='mjd').isot
    
        fits_keywords = {'DATE-OBS': date, **fits_keywords}
 
        return super().from_data(fits_keywords=fits_keywords, **columns)

    def update_uv(self):
        """

Update the (u, v) coordinates (UCOORD, VCOORD) using information of
the array and target information contained in OI_ARRAY and OI_TARGET 
tables.

Warnings
--------

(u, v) may differ from (UCOORD, VCOORD) determined by a data processing
software because

a. (UCOORD, VCOORD) can be averaged independently from MJD 
(see OIFITS standard)
b. Atmospheric refraction is dealt with approximately, while (UCOORD,
   VCOORD) may have none to full modelling of the atmosphere.


        """
        raise NotImplementedError('abstract class') 

    def get_sky_coord(self, max_distance=None):
        """
Get the full coordinates of the targets, including distance and 
motions, at the time (epoch) of observation.  Only the positions
are corrected for secular motion, the velocity and proper motions
are not computed.

Returns
-------

astropy.coordinates.SkyCoord object

        """
        # cross-reference rows of the Target HDU 
        thdu = self.get_targetHDU()
        indices = self._xmatch(refhdu=thdu, refname='TARGET_ID')
        refcoo =  thdu.get_sky_coord()

        coo = refcoo[indices]
        
        # apply epoch correction to positions
        obstime = _Time(self.MJD, format='mjd')
        coo = _coo.apply_space_motion(coo, obstime, correct_motion=False)

        if max_distance is not None:
            coo.distance[_np.isnan(coo.distance)] = max_distance

        return coo 

    def get_stauvw(self, *, refraction=False):
        """

Determine the (u, v, w) coordinates of each station with respect to
the centre of the array at the mean time of observation (MJD) using 
target and array information.

Arguments
---------

refraction (bool, optional, default: False)
    Whether the atmospheric refraction will be included in the
    calculation.

Returns
-------

uvw (float × NOBS × NSTA × 3)
    (u, v, w) coordinates for each observation and each station. NOBS
    is the number of observations (rows) in the table and NSTA the number
    of stations (NSTA = 3 for OI_T3, NSTA = 2 for OI_VIS and OI_VIS2,
    NSTA = 1 for OI_FLUX in uncalibrated mode)

Precision
---------

Sources of errors are:
1. About 7×10⁻⁵ per second of uncertainty on the effective mean time
of observation. 
2. About 5×10⁻⁵ at elevation of 45⁰ (goes as cot z), due to the 
uncertainty on temperature and pressure in the atmospheric refraction 
correction.
3. Uncertainty on the effective wavelength and relative humidity.

Discrepancies of a few centimetres for hectometric baselines have been 
observed for the VLTI.

Warning
-------

(u, v) may differ from the tabulated values in the tables (e.g. UCOORD, 
VCOORD) because data processing pipelines may average (UCOORD, VCOORD) 
independently from MJD (see OIFITS standard) and atmospheric refraction 
may be dealt with differently.

Raises
------

NotImplemented Error
    There is no associated OI_ARRAY table (possible in OIFITS v. 1)
        """   
        arrayHDU = self.get_arrayHDU() 
        waveHDU = self.get_wavelengthHDU()
        if arrayHDU is None:
            raise NotImplementedError('No OI_ARRAY is referenced')

        # If frame is SKY (e.g. pupil mask), STAXYZ are (u,v,w)
        # coordinates
        XYZ = self.get_staxyz()
        frame = arrayHDU.header.get('FRAME', 'GEOCENTRIC')
        if frame == 'SKY':
            return XYZ
 
        # enu = self.get_staenu()

        loc = arrayHDU.get_location()
        lat = loc.lat.to_value('rad')
        h = loc.height.value
        obswl = waveHDU.EFF_WAVE.mean() * _units.m
        waveHDU = self.get_wavelengthHDU()

        # Transform FK5 coordinates to apparent (atmosphere-refracted)
        # ITRS coordinates. We need to go through altaz because astropy.
        FK5 = self.get_sky_coord(max_distance=1 * _units.Mpc)
        altaz_frame = _coo.altaz_frame(loc, obswl=obswl, refraction=refraction)
        # print(f"{FK5.ra[0]=} {FK5.dec[0]=}")

        UVW = _np.empty_like(XYZ)      
 
        for i, (fk5, xyz) in enumerate(zip(FK5, XYZ)):

            altaz = fk5.transform_to(altaz_frame)
            altaz = _SkyCoord(altaz.replicate(pressure=0))
            itrs = altaz.itrs.spherical

            # (X, Y, Z)_ITRS -> (u, v, w)_SKY transform
            #
            # It is equivalent to the the (u, v, w) formula given by Eq. 2-30,
            # p. 25 in Synthesis Imaging in Radio Astronomy II, A Collection of 
            # Lectures from the Sixth NRAO/NMIMT Synthesis 
            # Imaging Summer School, G. B. Taylor, C. L. Carilli, and 
            # R. A. Perley. (eds.), ASP Conference Series, 180, 1999
            #
            # except they use a local frame centred on meridian and the 
            # hour angle, while here the frame is centred on the Greenwich
            # meridian and the longitude of the source in the ITRS is used. 
             
            lon = itrs.lon.value # target longitude in ITRS (opposite of hour
                                 # angle for an observer at Greenwich meridian)
            lat = itrs.lat.value # target declination in ITRS
            #if i == 0:
            #    print(f"{altaz.alt=} {altaz.az=}")
            #    print(f"{lon=} {lat=}") 
            wuv = _u.rotation3d(xyz, 'zy', [-lon, lat], degrees=True)
            uvw = _np.roll(wuv, -1, axis=-1)

            UVW[i] = uvw

        return UVW

# OIFITS1 Table 
class _DataHDU1(
        _DataHDU, 
        _MustHaveTargetHDU,
        _MayHaveArrayHDU,
        _MustHaveWavelengthHDU,
        _OIFITS1HDU
      ):
    _COLUMNS = [('TIME', True, '1D', (), None, 0., "s",
                    'Seconds since UT at start of observation')]
        

# OIFITS1 Table rev1
class _DataHDU11( 
         _DataHDU1,
         _OITableHDU1,
      ):

    def from_data(cls, *, fits_keywords={}, mjd, time=None, **columns):

        if time is None:
            time = (mjd - int(min(mjd))) * 24 * 3600 

        return super().from_data(fits_keywords=fits_keywords,
                                    mjd=mjd, time=time, **columns)

# OIFITS2 Table
class _DataHDU2(
        _DataHDU,
        _MustHaveTargetHDU,
        _MustHaveArrayHDU,
        _MustHaveWavelengthHDU,
        _MayHaveCorrHDU,
        _MayHaveInspolHDU,
        _OIFITS2HDU, 
      ):
    pass

# OIFITS2 Table rev1 (new table in OIFITS2)
class _DataHDU21(
        _DataHDU2,
        _OITableHDU1
      ):
    pass

# OIFITS2 Table rev2 (table was available in OIFITS1 but was updated)
class _DataHDU22(
        _DataHDU2,
        _OITableHDU2
      ):
    
    _COLUMNS = [('TIME', True, '1D', (), _u.is_zero, 0., "s",
                    'For backwards compatibility only')]

    @classmethod
    def from_data(cls, *, fits_keywords={}, **columns):

        columns['time'] = 0.

        return super().from_data(fits_keywords=fits_keywords, **columns)

