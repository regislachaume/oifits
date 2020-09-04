from .table import _OITableHDU, _OITableHDU11, _OITableHDU22
from .referenced import _Referenced
from astropy.coordinates import EarthLocation as _EarthLocation
from scipy.spatial.transform import Rotation as _rotation
from .. import utils as _u

import numpy as _np

class _ArrayHDUBase(_OITableHDU):

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
        if self.header['EXTNAME'] == 'OI_ARRAY':
            return self
        container = self.get_container()
        if container is None:
            return None
        return container.get_arrayHDU(arrname=self.get_arrname())

    def _verify(self, option='warn'):

        errors = super()._verify(option)
      
        # If no ARRNAME and STA_INDEX, fine (optional), but if
        # only one of them is present, there's a problem.
        sta_index = getattr(self, 'STA_INDEX', None)
        arrname = self.header.get('ARRNAME', None) 
        
        if not arrname:
            if sta_index is not None and not isinstance(self, _ArrayHDU):
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

    def __ge__(self, other):
        return (isinstance(self, _OITableHDU) and
                h1.get_arrname() == other.get('ARRNAME', ''))
    def __gt__(self, other):
        return self is not other and self >= other

class _MayHaveArrayHDU(_ArrayHDUBase):
    _CARDS = [('ARRNAME', False, _u.is_nonempty, None, 
        'Name of telescope array for cross-reference')]

class _MustHaveArrayHDU(_ArrayHDUBase):
    _CARDS = [('ARRNAME', True, _u.is_nonempty, None, 
        'Name of telescope array for cross-reference')] 

class _ArrayHDU(_MustHaveArrayHDU,_Referenced):
    
    _EXTNAME = 'OI_ARRAY'
    _REFERENCE_KEY = 'ARRNAME'
    _CARDS = [
        ('ARRAYX', True, _u.is_num, None, 'X coordinate of array centre'),
        ('ARRAYY', True, _u.is_num, None, 'Y coordinate of array centre'),
        ('ARRAYZ', True, _u.is_num, None, 'Z coordinate of array centre'),
    ]
    _COLUMNS = [
        ('TEL_NAME',  True, '<U16', (),   None,           None, None,
                'telescope name'), 
        ('STA_NAME',  True, '<U16', (),   None,           None, None,
                'station name'), 
        ('STA_INDEX', True, '>i2', (),   _u.is_strictpos, None, None,
                'station index for cross-reference'),
        ('DIAMETER',  True, '>f4', (),   _u.is_strictpos, None, "m",
                'effective diameter of the aperture'),
        ('STAXYZ',    True, '>f8', (3,), None,            None, "m",
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
        if len(_np.unique(sta_index)) == len(sta_index):
            return errors

        err_text = f"Repeated STA_INDEX in {type(self).__name__}"
        err = self.run_option(option, err_text, fixable=False)
        errors.append(err)

        return errors

    def __add__(self, other):

        return self.merge(other)

    def merge(self, *others):

        norm = _np.linalg.norm
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

    @classmethod
    def from_data(cls, arrname, *, version=2, 
            arrayxyz=None, lat=None, lon=None, alt=None, 
            frame='GEOCENTRIC', ellipsoid='WGS84',
            tel_name=None, sta_name=None, sta_index=None, diameter=None, 
            staxyz=None, sta_enu=None,
            fits_keywords={}, **columns):

        geodetic = lon is not None and lat is not None and alt is not None
        if geodetic:
            if frame == 'GEOCENTRIC':
                loc = _EarthLocation.from_geodetic(lon, lat, alt, ellipsoid)
                arrayxyz = [loc.x.value, loc.y.value, loc.z.value]
        if frame == 'SKY':
            arrayxyz = [0., 0., 0.]

        if lon is not None:
            fits_keywords['LON-OBS'] = (lon, 'longitude of array centre (deg)')
        if lat is not None:
            fits_keywords['LAT-OBS'] = (lat, 'latitude of array centre (deg)')
        if alt is not None:
            fits_keywords['ALT-OBS'] = (alt, 'altitude of array centre (m)')
   
        if sta_enu is not None:
            if frame == 'SKY':
                staxyz = sta_enu
            else:
                if not geodetic:
                    loc = _EarthLocation.from_geocentric(*arrayxyz, unit="m")
                    loc = loc.to_geodetic(ellipsoid=ellipsoid)
                    lon = loc.lon.value
                    lat = loc.lat.value
                rot = _rotation.from_euler('yz', [-lat, lon], degrees=True)
                sta_une = _np.roll(sta_enu, 1, axis=1)
                staxyz = rot.apply(sta_une)
            
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
    _CARDS = [('FRAME', True, _u._is_frame1, 'GEOCENTRIC', 
        'coordinate frame for array centre and stations')]

class ArrayHDU2(
        _ArrayHDU,
        _OITableHDU22, # OIFITS2, rev. 2
      ):
    _CARDS = [('FRAME', True, _u._is_frame2, 'GEOCENTRIC', 
        'coordinate frame for array centre and stations')]
    _COLUMNS = [
        ('FOV',     False, '>f8', (), _u.is_strictpos, None, "arcsec",
                'photometric field of view'), 
        ('FOVTYPE', False, '<U6', (), _u._is_fovtype2, None, None,
                'definition of the field of view')
    ]

new_array_hdu = _ArrayHDU.from_data
