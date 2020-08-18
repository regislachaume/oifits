from .base import _ValidHDU
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

class PrimaryHDU1(
        _PrimaryHDU
      ):
    _CARDS = [
        ('CONTENT', False, _u.is_oifits1, 'OIFITS1'),
    ]

class PrimaryHDU2(
        _PrimaryHDU
      ):
    _CARDS = [
        ('CONTENT',  True,  _u.is_oifits2,       'OIFITS2'),
        ('TELESCOP', True,  _u.is_nonempty,      'N/A'),
        ('INSTRUME', True,  _u.is_nonempty,      'N/A'),
        ('OBJECT',   True,  _u.is_nonempty,      'N/A'),
        ('REFERENC', False, _u.is_nonempty,      'N/A'),
        ('PROG_ID',  False, _u.is_nonempty,      'N/A'),
        ('PROCSOFT', False, _u.is_nonempty,      'N/A'),
        ('OBSTECH',  False, _u.is_nonempty,      'N/A'),
        ('RA',       False, _u.is_num,           None),
        ('DEC',      False, _u.is_num,           None),
        ('EQUINOX',  False, _u.is_num,           2000.0),
        ('RADECSYS', False, _u.is_nonempty,      'ICRS'),
        ('SPECSYS',  False, _u.is_nonempty,      'N/A'),
        ('TEXPTIME', False, _u.is_strictpos,     None),
        ('MJD-OBS',  False, _u.is_num,           None),
        ('MJD-END',  False, _u.is_num,           None),
        ('BASE_MIN', False, _u.is_strictposnum,  None),
        ('BASE_MAX', False, _u.is_strictposnum,  None),
        ('WAVELMIN', False, _u.is_strictposnum,  None),
        ('WAVELMAX', False, _u.is_strictposnum,  None),
        ('NUM_CHAN', False, _u.is_strictposint,  None),
        ('SPEC_RES', False, _u.is_strictposnum,  None),
        ('VIS2ERR',  False, _u.is_strictposnum,  None),
        ('VISPHERR', False, _u.is_strictposnum,  None),
        ('T3PHIERR', False, _u.is_strictposnum,  None),
    ]
