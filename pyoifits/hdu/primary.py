from .base import _ValidHDU, _OIFITS1HDU, _OIFITS2HDU
from .. import utils as _u

from astropy.io import fits as _fits

class _PrimaryHDU(
        _ValidHDU,
        _fits.PrimaryHDU
      ):
    def __repr__(self):
        return f"<{type(self).__name__} at {hex(id(self))} {self._diminfo()}>"
    def __str__(self):
        return f"<{type(self).__name__} {self._diminfo()}>"
    def _diminfo(self):
        if self.data is None:
            return "(void)"
        return f"({'×'.join(self.data.shape)})"

    @classmethod
    def match_header(cls, header):
        cards = getattr(cls, '_CARDS', None)
        if cards is None or 'CONTENT' not in cards['name']:
            return NotImplementedError
        version = cards[cards['name'] == 'CONTENT']['default'][0]
        return (_fits.PrimaryHDU.match_header(header) and
                header.get('CONTENT', 'OIFITS1') == version)

    @classmethod
    def __init_subclass__(cls):
        super().__init_subclass__()
        if hasattr(cls, '_CARDS') and 'CONTENT' in cls._CARDS['name']:
            _fits.hdu.base._BaseHDU.register_hdu(cls)

    def to_version(self, n):
        classes = {1: PrimaryHDU1, 2: PrimaryHDU2}
        cls = type(self)
        newcls = classes[n]
        if cls == newcls:
            return self
        hdu = newcls(data=self.data, header=self.header)
        hdu.header['CONTENT'] = f"oifits{n}"
        return hdu
    
    def __init__(self, data=None, header=_fits.Header(), keywords={}):
        
        keywords = {k.upper(): v for k, v in keywords.items() if v is not None}

        if 'content' in keywords:
            del keywords['content']
    
        for card in self._CARDS:
            name = card['name']
            comment = card['comment']
            if name in keywords:
                value = keywords[name]
                if isinstance(value, (list, tuple)):
                    value, comment = value[0:2]
                header.set(name, value, comment)
                del keywords[name]
            elif card['required'] and name not in header:
                value = card['default']
                if value is not None:
                    header.set(name, value, comment)
        
        for name, value in keywords.items():
            if not isinstance(value, (tuple, list)):
                value = [value]
            header.set(name, *value)
        
        super().__init__(data=data, header=header)
                


class PrimaryHDU1(
        _PrimaryHDU,
        _OIFITS1HDU,
      ):
    _CARDS = [
        ('CONTENT', False, _u.is_oifits1, 'OIFITS1', 'format by Pauls et al. (2005), PASP 117,1125'),
    ]

        

class PrimaryHDU2(
        _PrimaryHDU,
        _OIFITS2HDU,
      ):
    _CARDS = [
        ('CONTENT',  True,  _u.is_oifits2,       'OIFITS2', 'format by Duvert et al. (2017), A&A 597, A8'),
        ('TELESCOP', True,  _u.is_nonempty,      'N/A', 'name of the telescope array'),
        ('INSTRUME', True,  _u.is_nonempty,      'N/A', 'name of the interferometric instrument'),
        ('OBJECT',   True,  _u.is_nonempty,      'N/A', 'astronomical object'),
        ('REFERENC', False, _u.is_nonempty,      'N/A', 'bibliographic reference'),
        ('PROG_ID',  False, _u.is_nonempty,      'N/A', 'observing programme ID'),
        ('PROCSOFT', False, _u.is_nonempty,      'N/A', 'data processing software'),
        ('OBSTECH',  False, _u.is_nonempty,      'N/A', 'beam recombination technique'),
        ('RA',       False, _u.is_num,           None, 'right ascension (deg)'),
        ('DEC',      False, _u.is_num,           None, 'declination (deg)'),
        ('EQUINOX',  False, _u.is_num,           2000., 'equinox (yr)'),
        ('RADECSYS', False, _u.is_nonempty,      'ICRS', 'celestial coordinate frame'),
        ('SPECSYS',  False, _u.is_nonempty,      'N/A', ''),
        ('TEXPTIME', False, _u.is_strictpos,     None, 'total exposure time'),
        ('MJD-OBS',  False, _u.is_num,           None, 'MJD at start of observation (d)'),
        ('MJD-END',  False, _u.is_num,           None, 'MJD at end of observation (d)'),
        ('BASE_MIN', False, _u.is_strictposnum,  None, 'minimum baseline length (m)'),
        ('BASE_MAX', False, _u.is_strictposnum,  None, 'maximum baseline length (m)'),
        ('WAVELMIN', False, _u.is_strictposnum,  None, 'minimum wavelength (m)'),
        ('WAVELMAX', False, _u.is_strictposnum,  None, 'maximum wavelength (m)'),
        ('NUM_CHAN', False, _u.is_strictposint,  None, 'number of spectral channels'),
        ('SPEC_RES', False, _u.is_strictposnum,  None, 'spectral resolution'),
        ('VIS2ERR',  False, _u.is_strictposnum,  None, 'typical uncertainty in squared visibilities'),
        ('VISPHERR', False, _u.is_strictposnum,  None, 'typical uncertainty in phases (deg)'),
        ('T3PHIERR', False, _u.is_strictposnum,  None, 'typical uncertainty in closure phases (deg)'),
    ]
