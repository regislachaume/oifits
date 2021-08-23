"""Implementation of the OI_ARRAY binary table extension."""
import numpy as np
from astropy.coordinates import EarthLocation 

from .. import utils as u
from .table import _OITableHDU, _OITableHDU11, _OITableHDU22
from .referenced import _Referenced

__all__ = ["ArrayHDU1", "ArrayHDU2", "new_array_hdu"] 

class _ArrayHDUBase(_OITableHDU):

    def _get_array_field(self, name, shape='none', flatten=False,
            concatenate=False, default=None):

        refhdu = self.get_arrayHDU()

        if refhdu is None:
            if default is None:
                return None
            val = default
        else:
            val = self._xmatch(refhdu, 'STA_INDEX', 
                            name=name, concatenate=concatenate)

        return self._resize_data(val, shape, flatten)

    def get_sta_name(self, shape='none', flatten=False):
        return self._get_array_field('STA_NAME', shape, flatten)

    def get_staxyz(self):
        return self._get_array_field('STAXYZ')

    def get_staenu(self):
        """

Determine the (East, North, up) coordinates of each station with respect
with the centre of the array

Arguments
---------

None

Returns
-------

uvw (float × NOBS × NSTA × 3)
    (u, v, w) coordinates for each observation and each station. NOBS
    is the number of observations (rows) in the table and NSTA the number
    of stations (NSTA = 3 for OI_T3, NSTA = 2 for OI_VIS and OI_VIS2,
    NSTA = 1 for OI_FLUX in uncalibrated mode)

Raises
------

NotImplementedError 
    The array is not defined in the GEOCENTRIC frame

        """
        xyz = self.get_staxyz()
        frame = self.header.get('FRAME', 'GEOCENTRIC')

        if frame != 'GEOCENTRIC':
            txt = f"(E, N, U) coordinates not implemented for frame {frame}"
            raise NotImplementedErrror(txt)

        # (X, Y, Z)_ITRS -> (E, N, U)_local 
        loc = self.get_location()
        lon = loc.lon.to_value('rad')
        lat = loc.lat.to_value('rad')
        staune = u.rotation3d(xyz, 'zy', [-lon, lat], degrees=False)
        staenu = np.roll(staune, -1, axis=-1)
       
        return staenu
 
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

    def get_location(self):
        
        arrayHDU = self.get_arrayHDU()
        header = arrayHDU.header
        xyz = [header.get('ARRAY' + var, 0) for var in ['X', 'Y', 'Z']]
        loc = EarthLocation.from_geocentric(*xyz, unit="m")

        return loc

    def get_arrayHDU(self):
        if self.header['EXTNAME'] == 'OI_ARRAY':
            return self
        container = self.get_container()
        if container is None:
            return None
        return container.get_arrayHDU(arrname=self.get_arrname())

    def _verify(self, option='warn'):
        
        errors = super()._verify(option)
        
        sta_index = getattr(self, 'STA_INDEX', None)
        arrname = self.header.get('ARRNAME', None) 
        
        # If no ARRNAME and STA_INDEX, fine when optional, but if
        # only one of them is present, there's a problem.

        loc = f"{type(self).__name__} object"
        
        if arrname is None and sta_index is not None:

            text = f"STA_INDEX column present with no ARRNAME in {loc}"
            err = self.run_option(option, text, fixable=False)
            errors.append(err)
         
        elif arrname is not None and sta_index is None:

            text = f"ARRNAME is given but not STA_INDEX column in {loc}"
            err = self.run_option(option, text, fixable=False)
            errors.append(err)

        return errors

    def __ge__(self, other):
        return (isinstance(self, _OITableHDU) and
                h1.get_arrname() == other.get('ARRNAME', ''))
    def __gt__(self, other):
        return self is not other and self >= other


class _MayHaveArrayHDU(_ArrayHDUBase):
    _CARDS = [('ARRNAME', False, u.is_nonempty, None, 
        'Name of telescope array for cross-reference')]

class _MustHaveArrayHDU(_ArrayHDUBase):
    _CARDS = [('ARRNAME', True, u.is_nonempty, None, 
        'Name of telescope array for cross-reference')] 

class _ArrayHDU(_MustHaveArrayHDU,_Referenced):
    
    _EXTNAME = 'OI_ARRAY'
    _REFERENCE_KEY = 'ARRNAME'
    _CARDS = [
        ('ARRAYX', True, u.is_num, None, 'X coordinate of array centre'),
        ('ARRAYY', True, u.is_num, None, 'Y coordinate of array centre'),
        ('ARRAYZ', True, u.is_num, None, 'Z coordinate of array centre'),
    ]
    _COLUMNS = [
        ('TEL_NAME',  True, '16A', (),  None,            None, None,
                'telescope name'), 
        ('STA_NAME',  True, '16A', (),  None,            None, None,
                'station name'), 
        ('STA_INDEX', True, '1I',  (),  u.is_int,       None, None,
                'station index for cross-reference'),
        ('DIAMETER',  True, '1E', (),   u.is_strictpos, None, "m",
                'effective diameter of the aperture'),
        ('STAXYZ',    True, '3D', (3,), None,            None, "m",
                'station coordinates relative to array centre'),
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
        if len(np.unique(sta_index)) == len(sta_index):
            return errors

        err_text = f"Repeated STA_INDEX in {type(self).__name__}"
        err = self.run_option(option, err_text, fixable=False)
        errors.append(err)

        return errors

    def __add__(self, other):

        return self.merge(other)

    def merge(self, *others):

        norm = np.linalg.norm
        dist_max = self.get_container()._merge_station_distance

        def eq(x, y):
            return ((x['STA_NAME'] == y['STA_NAME']) and
                    (x['TEL_NAME'] == y['TEL_NAME']) and
                    (norm(x['STAXYZ'] - y['STAXYZ']) <= dist_max))
        
        return self._merge_helper(*others, id_name='STA_INDEX', equality=eq)

    def __mod__(self, other):

        h1, h2 = self.header, other.header
        dist_max = self.get_container()._merge_array_distance
        return (self & other and
                h1['ARRNAME'] == h2['ARRNAME'] and
                h1['FRAME'] == h2['FRAME'] and
                abs(h1['ARRAYX'] - h2['ARRAYX']) <= dist_max and
                abs(h1['ARRAYY'] - h2['ARRAYY']) <= dist_max and
                abs(h1['ARRAYZ'] - h2['ARRAYZ']) <= dist_max)

    def is_referred_to_by(self, other):
        return (not isinstance(other, _ArrayHDU) and
                isinstance(other, _ArrayHDUBase) and
                self.get_arrname() == other.get_arrname())
    
    def reindex(self):

        # cannot reindex a standalone table
        if self.get_container() is None:
            raise RuntimeError('cannot reindex a standalone OI_ARRAY')

        hdus = [self, *self.get_referrers()]

        for new, old in enumerate(self.STA_INDEX, start=1):
            for h in hdus:
                h.STA_INDEX[h.STA_INDEX == old] = new

    @classmethod
    def from_data(cls, *, arrname, version=2, 
            arrayxyz=None, lat=None, lon=None, alt=None, 
            frame='GEOCENTRIC', ellipsoid='WGS84',
            tel_name=None, sta_name=None, sta_index=None, diameter=None, 
            staxyz=None, staenu=None,
            fits_keywords={}, **columns):
        """

Build an OI_ARRAY extension from data. In the following, NSTA is the
number of stations.

Arguments
---------

version:
    OIFITS version if it cannot be determined from context (optional)
ellipsoid (str)
    Geodetic coordinate system (optional, default is WGS84)
calibrated (bool)
    Whether flux are calibrated or uncalibrated (measured at each station)

arrname (str)
    name of the array
frame (str)
    frame (optional, default is GEOCENTRIC)
fits_keywords (dict)
    additional FITS keywords (optional)
arrayxyz (float, 3) 
or lat, lon, alt (float):
    Coordinates of the centre of the array. 
    * arrayxyz uses the geocentric frame in metres.
    * lat, lon, alt are the geodetic coordinates in degrees and metres.
    Only for GEOCENTRIC frame.

sta_index (int, NSTA)
    Station index for cross-reference
tel_name (str, NSTA)
    Names of the telescopes
sta_name (str, NSTA)
    Names of the stations
diameter (float, NSTA)
    Effective aperture of the telescopes
staxyz (float, NSTA × 3)
or staenu (float, NSTA × 3)
    Coordinates of the stations relative to array centre. staenu uses
    the local frame (East-North-up) and sta_xyz the specified frame
    (Cartesian geocentric if GEOCENTRIC, local frame if SKY). In metres.

Arguments in OIFITS2 only
-------------------------

fov (float)
    field of view in arcsec
fovtype (str)
    field of view type is FWHM or RADIUS 

Additional arguments
--------------------

Any additional keyword argument will be appended as a non-standard FITS 
column with its name prefixed with NS_ 

        """
        geodetic = lon is not None and lat is not None and alt is not None
        if geodetic:
            if frame == 'GEOCENTRIC':
                loc = EarthLocation.from_geodetic(lon, lat, alt, ellipsoid)
                arrayxyz = [loc.x.value, loc.y.value, loc.z.value]
        if frame == 'SKY':
            arrayxyz = [0., 0., 0.]

        if lon is not None:
            fits_keywords['LON-OBS'] = (lon, 'longitude of array centre (deg)')
        if lat is not None:
            fits_keywords['LAT-OBS'] = (lat, 'latitude of array centre (deg)')
        if alt is not None:
            fits_keywords['ALT-OBS'] = (alt, 'altitude of array centre (m)')
   
        if staenu is not None:
            if frame == 'SKY':
                staxyz = staenu
            else:
                if not geodetic:
                    loc = EarthLocation.from_geocentric(*arrayxyz, unit="m")
                    loc = loc.to_geodetic(ellipsoid=ellipsoid)
                    lon = loc.lon.value
                    lat = loc.lat.value
                # (E, N, U) -> (X, Y, Z)_ITRS transform
                staune = np.roll(staenu, 1, axis=-1)
                staxyz = u.rotation3d(staune, 'yz', [-lat, lon], degrees=True)
        if sta_index is None:
            sta_index = list(range(1, len(tel_name) + 1))
        
        fits_keywords = dict(arrname=arrname, frame=frame,
            arrayx=arrayxyz[0], arrayy=arrayxyz[1], arrayz=arrayxyz[2],
            **fits_keywords)
   
        columns = dict(tel_name=tel_name, sta_name=sta_name, 
            sta_index=sta_index, diameter=diameter, staxyz=staxyz,
            **columns)
        
        return super().from_data(version=version, 
                    fits_keywords=fits_keywords, **columns)

class ArrayHDU1(
        _ArrayHDU,
        _OITableHDU11, # OFITS1, rev. 1
      ):
    """

First revision of the OI_ARRAY binary table, OIFITS v. 1

    """
    _CARDS = [('FRAME', True, u.is_frame1,  'GEOCENTRIC', 
        'coordinate frame for array centre and stations')]

class ArrayHDU2(
        _ArrayHDU,
        _OITableHDU22, # OIFITS2, rev. 2
      ):
    """

Second revision of the OI_ARRAY binary table, OIFITS v. 2

    """
    _CARDS = [('FRAME', True, u.is_frame2,  'GEOCENTRIC', 
        'coordinate frame for array centre and stations')]
    _COLUMNS = [
        ('FOV',     False, '1D', (), u.is_strictpos, None, "arcsec",
                'photometric field of view'), 
        ('FOVTYPE', False, '6A', (), u.is_fovtype2,  None, None,
                'definition of the field of view')
    ]
        
    def _verify(self, option='warn'):

        errors = super()._verify(option)

        # Verify STA_INDEX is correct (>= 1)

        indices = np.unique(self.STA_INDEX)
        if not all(indices >= 1):

            err_text = "'STA_INDEX' should be ≥ 1."
            fix_text = "Substituting with values >= 1"
            fix = lambda h=self: self.reindex()
            
            if self.get_container() is not None:
                err = self.run_option(option, err_text, fix_text, fix=fix)
            else:
                err = self.run_option(option, err_text, fixable=False)
            
            errors.append(err)
        
        return errors

new_array_hdu = _ArrayHDU.from_data
