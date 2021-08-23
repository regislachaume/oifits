from astropy.io import fits as fits
import numpy as np

from .base import _ValidHDU, _OIFITS1HDU, _OIFITS2HDU
from .. import utils as u
       
__all__ = ["PrimaryHDU1", "PrimaryHDU2", "new_primary_hdu"]
     
# mean of medians, where X = [x1, ..., xn]
# testing > 0 also exclude NULL values (NaNs)

def _mean_med(X):
    X = [x_ok for x in X if len(x_ok := x[x>0])]
    if X:
        return np.mean([np.ma.median(x) for x in X])


class _PrimaryHDU(
        _ValidHDU,
        fits.PrimaryHDU
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
        return (fits.PrimaryHDU.match_header(header) and
                header.get('CONTENT', 'OIFITS1') == version)

    @classmethod
    def __init_subclass__(cls):
        super().__init_subclass__()
        if hasattr(cls, '_CARDS') and 'CONTENT' in cls._CARDS['name']:
            fits.hdu.base._BaseHDU.register_hdu(cls)

    def _to_version(self, n):

        cls = type(self)

        classes = cls.__base__.__subclasses__()
        for newcls in classes:
            if newcls._OI_VER == n:
                break
        else:
            raise ValueError(f'unknown OIFITS version number: {n}')

        primary = self.copy()

        if cls == newcls:
            return primary

        hdu = newcls(data=primary.data, header=primary.header)
        hdu.header['CONTENT'] = f"OIFITS{n}"

        return hdu
    
    def __init__(self, data=None, header=None, keywords={}):
      
        if header is None:
            header = fits.Header()

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

    def _update_header(self):

        from astropy.time import Time

        cont = self.get_container()
        if not cont:
            return
        
        header = self.header
        targets = cont.get_targetHDU().data
        datahdus = cont.get_dataHDUs()

        # ensure content is just after the EXTEND keyword
        cards = self._CARDS 
        card = cards[cards['name'] == 'CONTENT'][0]
        if 'CONTENT' in header:
            del header['CONTENT']
        value = ('CONTENT', card['default'], card['comment'])
        header.insert('EXTEND', value, after=True)

        # fix dates
        mjdobs = min(h.MJD.min() for h in datahdus)
        mjdend = max(h.MJD.max() for h in datahdus)
        header['DATE-OBS'] = Time(mjdobs, format='mjd').isot[0:19]
        header['DATE'] = Time.now().isot[0:19]

        # single target    
        if len(targets) == 1:
    
            header['OBJECT'] = targets['TARGET'][0]
            header['RA'] = targets['RAEP0'][0]
            header['DEC'] = targets['DECEP0'][0]
            header['EQUINOX'] = targets['EQUINOX'][0]
            header['MJD-OBS'] = mjdobs
            header['MJD-END'] = mjdend

        else:

            for keyw in ['RA', 'DEC', 'UTC', 'LST', 'EQUINOX', 'RADECSYS',
                'MJD-OBS', 'MJD-END']:
                if keyw in header:
                    del header[keyw]

    def update_header(self):
        
        self._update_header()
        self._fill_missing_keys()

    def _fill_missing_keys(self):

        # missing keyword
        cards = self._CARDS
        header = self.header

        for card in cards[cards['required']]:
            name, default = card['name'], card['default']
            if name not in header and default is not None:
                header[name] = default

        # update header comments
        for name, comment in zip(cards['name'], cards['comment']):
            if name in header:
                header.set(name, comment=comment)


class PrimaryHDU1(
        _PrimaryHDU,
        _OIFITS1HDU,
      ):

    _CARDS = [
        ('CONTENT',  False, u.is_oifits1,   'OIFITS1', 'format by Pauls et al. (2005), PASP 117,1125'),
        ('TELESCOP', False, u.is_nonempty,     'N/A',  'name of the telescope array'),
        ('INSTRUME', False, u.is_nonempty,     'N/A',  'name of the interferometric instrument'),
        ('ORIGIN',   False, u.is_nonempty,     'N/A',  'institution responsible of file creation'),
        ('OBSERVER', False, u.is_nonempty,     'N/A',  'person who took the data'),
        ('OBJECT',   False,  u.is_nonempty,    'N/A',  'target designation'),
        ('DATE',     False, u.is_date,         'N/A',  'date the HDU was written'),
        ('DATE-OBS', False, u.is_date,         'N/A',  'date of start of the observations'),
        ('RA',       False, u.is_num,           None,  'right ascension (deg)'),
        ('DEC',      False, u.is_num,           None,  'declination (deg)'),
        ('EQUINOX',  False, u.is_num,           2000., 'equinox (yr)'),
        ('RADECSYS', False, u.is_nonempty,     'ICRS', 'celestial coordinate frame'),
        ('MJD-OBS',  False, u.is_num,           None, 'MJD at start of observation (d)'),
        ('MJD-END',  False, u.is_num,           None, 'MJD at end of observation (d)'),
    ]


class PrimaryHDU2(
        _PrimaryHDU,
        _OIFITS2HDU,
      ):

    _CARDS = [
        ('CONTENT',  True,  u.is_oifits2,       'OIFITS2', 'format by Duvert et al. (2017), A&A 597, A8'),
        ('TELESCOP', True,  u.is_nonempty,      'N/A', 'name of the telescope array'),
        ('INSTRUME', True,  u.is_nonempty,      'N/A', 'name of the interferometric instrument'),
        ('ORIGIN',   True,  u.is_nonempty,      'N/A', 'institution responsible of file creation'),
        ('OBSERVER', True,  u.is_nonempty,      'N/A', 'person who took the data'),
        ('OBJECT',   True,  u.is_nonempty,      'N/A', 'astronomical object'),
        ('DATE',     True,  u.is_date,          'N/A', 'date the HDU was written'),
        ('DATE-OBS', True,  u.is_date,          'N/A', 'date of start of the observations'),
        ('RA',       False, u.is_num,           None,  'right ascension (deg)'),
        ('DEC',      False, u.is_num,           None,  'declination (deg)'),
        ('EQUINOX',  False, u.is_num,           2000., 'equinox (yr)'),
        ('RADECSYS', False, u.is_nonempty,     'ICRS', 'celestial coordinate frame'),
        ('MJD-OBS',  False, u.is_num,           None, 'MJD at start of observation (d)'),
        ('MJD-END',  False, u.is_num,           None, 'MJD at end of observation (d)'),

        ('INSMODE',  True,  u.is_nonempty,      'N/A', 'instrument mode'),
        ('REFERENC', False, u.is_nonempty,      'N/A', 'bibliographic reference'),
        ('PROG_ID',  False, u.is_nonempty,      'N/A', 'observing programme ID'),
        ('PROCSOFT', False, u.is_nonempty,      'N/A', 'data processing software'),
        ('OBSTECH',  False, u.is_nonempty,      'N/A', 'beam recombination technique'),
        ('SPECSYS',  False, u.is_nonempty,      None,  'reference frame for spectral coord.'),
        ('TEXPTIME', False, u.is_strictpos,     None, 'ellapsed time during observation (s)'),
        ('BASE_MIN', False, u.is_strictposnum,  None, 'minimum baseline length (m)'),
        ('BASE_MAX', False, u.is_strictposnum,  None, 'maximum baseline length (m)'),
        ('WAVELMIN', False, u.is_strictposnum,  None, 'minimum wavelength (m)'),
        ('WAVELMAX', False, u.is_strictposnum,  None, 'maximum wavelength (m)'),
        ('NUM_CHAN', False, u.is_strictposint,  None, 'number of spectral channels'),
        ('SPEC_RES', False, u.is_strictposnum,  None, 'spectral resolution'),
        ('VIS2ERR',  False, u.is_strictposnum,  None, 'typ. uncertainty in squared vis. amp. (%)'),
        ('VISPHERR', False, u.is_strictposnum,  None, 'typ. uncertainty in phases (deg)'),
        ('T3PHIERR', False, u.is_strictposnum,  None, 'typ. uncertainty in closure phases (deg)'),
    ]

    def _update_header(self):

        import re

        super()._update_header()
        
        cont = self.get_container()
        if not cont:
            return

        header = self.header
        targets = cont.get_targetHDU().data
        datahdus = cont.get_dataHDUs()
        wavehdus = cont.get_wavelengthHDUs()

        if len(targets) == 1:

            vis2hdus = cont.get_vis2HDUs()
            vishdus = cont.get_visHDUs()
            t3hdus = cont.get_t3HDUs()

            wavelmin = min(h.EFF_WAVE.min() for h in wavehdus)
            header['WAVELMIN'] = float(f"{wavelmin:10.4g}")
            wavelmax = max(h.EFF_WAVE.max() for h in wavehdus)
            header['WAVELMAX'] = float(f"{wavelmax:10.4g}")

            if len(wavehdus) == 1:

                w = wavehdus[0]
                header['NUM_CHAN'] = len(w.data)
                res = np.mean(w.EFF_WAVE / w.EFF_BAND)
                header['SPEC_RES'] = round(res, 1)

            # baselines

            uv2 = [d.get_uv() ** 2 for d in datahdus]
            b = np.hstack([np.sqrt(x[::2] + x[1::2]).ravel() for x in uv2])
            header['BASE_MIN'] = round(b.min(), 2)
            header['BASE_MAX'] = round(b.max(), 2)

            # typical uncertainties (taking mean of HDU medians where bad
            # errors < 0 or NaNs are removed)
            
            eps = 1e-10 # avoid division by 0
            err = _mean_med([h.VIS2ERR / abs(eps + abs(h.VIS2DATA)) 
                                                        for h in vis2hdus])
            if err:
                header['VIS2ERR'] = float(f"{err:10.4f}")

            err = _mean_med([h.T3PHIERR for h in t3hdus])
            if err:
                header['T3PHIERR'] = float(f"{err:10.4f}")
            
            err = _mean_med([h.VISPHIERR for h in vishdus])
            if err:
                header['VISPHERR'] = float(f"{err:10.4f}")

            # TEXPTIME. Should be the same, but because some interferograms
            # get tossed for one observable and not the other, there may
            # be some variations in INT_TIME.  Taking the max, per standard. 
            t = np.unique([t for hdu in datahdus for t in hdu.INT_TIME])
            if t[-1] - t[0] < 0.3 * t[0]:
                header['TEXPTIME'] = t[-1]

        else:

            header['OBJECT'] = 'MULTI'
            
            for keyw in ['WAVELMIN', 'WAVELMAX', 'NUM_CHAN', 'SPEC_RES',
                'TEXPTIME', 'BASE_MIN', 'BASE_MAX',
                'VIS2ERR', 'VISPHERR', 'T3PHIERR']:
                if keyw in header:
                    del header[keyw]

        arrnames = np.unique([n for h in datahdus if (n := h.get_arrname())])
        header['TELESCOP'] = arrnames[0] if len(arrnames) == 1 else 'MULTI'
        
        insnames = np.unique([h.get_insname() for h in datahdus])
       
        # INSTRUME is harder to guess, as there is no norm about how INSNAME is
        # given.   
        ins = [re.sub('([A-Za-z]+).*', '\\1', i).upper() for i in insnames]
        ins = np.unique(ins)
        header['INSTRUME'] = ins[0] if len(ins) == 1 else 'MULTI'

        # INSMODE 
        header['INSMODE'] = insnames[0] if len(insnames) == 1 else 'MULTI'

def new_primary_hdu(version=2, *, header=None, data=None):
    cls = _PrimaryHDU.get_class(version=version) 
    if header is not None:
        header = header.copy()
    return cls(data=None, header=None)
    
