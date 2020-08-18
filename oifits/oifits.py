from astropy.io import fits as _fits
from astropy import table as _table
from astropy.time import Time as _Time
from numpy import ma as _ma
import numpy as _np
import re as _re 

from .hdu.base import _ValidHDU
from .hdu.table import _OITableHDU
from .hdu.data import _DataHDU
from .hdu.target import _TargetHDU, TargetHDU1, TargetHDU2
from .hdu.array import _ArrayHDU, ArrayHDU1, ArrayHDU2
from .hdu.wavelength import _WavelengthHDU, WavelengthHDU1, WavelengthHDU2
from .hdu.corr import _CorrHDU, CorrHDU1
from .hdu.inspol import _InspolHDU, InspolHDU1
from .hdu.t3 import T3HDU1, T3HDU2
from .hdu.vis import VisHDU1, VisHDU2
from .hdu.vis2 import Vis2HDU1, Vis2HDU2
from .hdu.flux import FluxHDU1
from .hdu.primary import _PrimaryHDU, PrimaryHDU1, PrimaryHDU2
from . import utils as _u


def open(filename, mode='readonly', lazy_load_hdus=True, **kwargs):
    # we need to read the primary HDU to see which FITS version it is
    if mode == 'ostream':
        cls = OIFITS2
    else:
        with _fits.open(filename, lazy_load_hdus=True) as hdulist:
            if list.__len__(hdulist): # len(h) loads all HDUs
                content = hdulist[0].header.get('CONTENT', '')
                cls = OIFITS2 if content == 'OIFITS2' else OIFITS1
            else:
                cls = OIFITS2
    hdus = cls.fromfile(filename, lazy_load_hdus=lazy_load_hdus, **kwargs)
    return hdus

def openlist(filenames):

    hdulists = [open(f, lazy_load_hdus=False) for f in filenames]
    hdulist = _merge(*hdulists, _inplace=True)

    return hdulist

def merge(*hdulists):

   return _merge(*hdulists, _inplace=False)     

def _merge(*hdulists, _inplace=False):

    # We need to copy everything once, because merging will need to change
    # a lot of stuff in the tables (TARGETID, STA_INDEX, ARRNAME, etc.)
    # It's inefficient if some big tables can be merged after, but no
    # easy cop out.
    if not _inplace:
        hdulists = [hdulist.copy() for hdulist in hdulists]

    for hdulist in hdulists:
        hdulist.verify('silentfix+ignore')

    # Use the latest OIFITS version used by the OIFITS to be merged
    oirev = [getattr(hdulist, '_OI_VER', 0) for hdulist in hdulists]
    cls = type(hdulists[_np.argmax(oirev)])
    
    # A flat array containing all HDUs.  They keep knowledge of their
    # container.
    hdus = [hdu for hdulist in hdulists for hdu in hdulist]

    # The processing steps:
    # * build a composite header from OIFITS primary headers
    # * convert non-void primary headers into image HDUs
    _process_primaryHDUs(hdus)

    # * remove duplicates (OI tables only)
    _remove_equal_OITableHDUs(hdus)

    # * merge TargetHDU, ArrayHDU when possible
    _merge_OITableHDUs(hdus, cls=(_TargetHDU, _ArrayHDU))

    # * duplicate references ARRNAME, INSNAME, CORRNAME leads
    #   to renaming.  Unfortunately HDUs to be merged will be
    #   needlessly copied here
    _rename_conflicting_OITableHDUs(hdus)

    # * merge interferometric data tables where possible
    _merge_OITableHDUs(hdus, cls=_DataHDU)

    # * update the primary HDU header inferring some keywords from 
    #   contents
    merged = cls(hdus)

    merged.update_primary_header()
    
    # * final tweaks
    merged.sort()
    merged.update_extver()

    return merged

def _process_primaryHDUs(hdus):
    
    # merge primary headers of OIFITs primary HDUs 
    headers = [hdu.header for hdu in hdus if isinstance(hdu, _PrimaryHDU)]
    header = _u.merge_fits_headers(*headers)
    hdus[0].header = header
    # delete void OIFIT Primary headers and convert others to images
    for i in range(len(hdus)-1, 0, -1):
        hdu = hdus[i]
        if isinstance(hdu, _fits.PrimaryHDU):
            if hdu.data is not None:
                hdus[i] = _fits.ImageHDU(hdu)
            else:
                del hdus[i]
       
def _remove_equal_OITableHDUs(hdus, cls=_OITableHDU):
    
    for i, hdu in enumerate(hdus):

        if not isinstance(hdu, cls):
            continue

        for j in range(len(hdus)-1, i, -1):
            if hdus[j] == hdu:
                del hdus[j]     


def _rename_conflicting_OITableHDUs(hdus, 
        cls=(_WavelengthHDU, _ArrayHDU, _CorrHDU)):

    for i, hdu1 in enumerate(hdus):

        if not isinstance(hdu1, cls):
            continue
            
        refkey = getattr(hdu1, '_REFERENCE_KEY', None)
        if not refkey:
            continue 

        # We then look if any further HDU has a duplicate name of hdus[i]. 
        # Sequential Name behaves like a string 'refname_number' aware
        # of its own numbering.

        hdus2 = hdus[i + 1:]
        
        refname1 = _u.SequentialName(hdu1.header[refkey])
        refhdus2 = [hdu2 for hdu2 in hdus2 if hdu2 & hdu1]
        refnames2 = [_u.SequentialName(h.header[refkey]) for h in refhdus2] 
        equal = [refname1 == refname2 for refname2 in refnames2]
        indices = _np.argwhere(equal)[:,0]
        nequal = len(indices)
        if not nequal:
            continue

        # If there are equal names, then rename using name_number.  We
        # need to look at name_number to pick a non used number...
        # seqname1 % seqname2 means same name but different number.  From
        # the list, we can pick available names. 
     
        refs = [refname2 for refname2 in refnames2 if refname2 % refname1] 
        new_refnames2 = refname1.next_available(refs, nequal)

        for index in indices:

            refname2 = refnames2[index]
            refhdu2 = refhdus2[index]
            new_refname2 = str(new_refnames2.pop(0))

            container = refhdu2.get_container()
            if not container:
                continue
        
            referrers = container.get_referrers(refhdu2)
            for h in referrers:
                h.header[refkey] = new_refname2
 
            refhdu2.header[refkey] = new_refname2
            
            if refkey == 'INSNAME':
                for h in container.getInspolHDU():
                    insname = h.data['INSNAME'] 
                    h.data['INSNAME'][insname == refname2] = new_refname2

def _merge_OITableHDUs(hdus, cls=_OITableHDU):

    # find all equal HDUs and remove them

    for i, hdu in enumerate(hdus):    

        if not isinstance(hdu, cls):
            continue
        
        # Find all mergeable HDUs, merge them 
        nhdu = len(hdus)
        is_mergeable = [hdus[k] % hdu if k > i else False for k in range(nhdu)]
        where_mergeable = _np.argwhere(is_mergeable)[:,0]
        if len(where_mergeable):
            mergeable = [hdus[m] for m in where_mergeable]
            hdus[i] = hdus[i].merge(*mergeable)
            for j in where_mergeable[::-1]:
                del hdus[j]
    
def set_merge_settings(*, station_distance=0.1, array_distance=10, 
       target_distance=2.5e-8, target_name_match=False):
        _OIFITS._merge_station_distance = station_distance
        _OIFITS._merge_array_distance = array_distance
        _OIFITS._merge_target_distance = target_distance
        _OIFITS._merge_target_name_match = target_name_match

class _OIFITS(_fits.HDUList):
    """Top-level class of Optical Interferometry FITS format.
When a fits file is opened a HDUList object is returned."""

    _merge_station_distance = 0.1    # 0.1 m
    _merge_array_distance = 10       # 10 m
    _merge_target_distance = 2.5e-8  # ~0.005 mas 
    _merge_target_name_match = False # allow different target designations.


    def __repr__(self):

        name = type(self).__name__
        str_ = [str(h) for h  in self]
        return f"<{name} at {hex(id(self))}: {' '.join(str_)}>"

    def __str__(self):

        name = type(self).__name__
        str_ = [str(h) for h  in self]

        return f"<{name}: {primary} {' '.join(str_)}>"
   
    def __init__(self, hdus=[], file=None, *, _copy_hdus=True):
    
        super().__init__(hdus=hdus, file=file)  
        for hdu in list.__iter__(hdus):
            hdu._container = self

    # original _read_next_hdu() uses super().append(), ruining any clean 
    # attempt to subclass HDUList
    def _read_next_hdu(self):
        
        has_new_hdu = super()._read_next_hdu()
        if has_new_hdu:
            last_index = list.__len__(self) - 1 # len(x) would load all HDUs
            hdu = self[last_index]
            if isinstance(hdu, _OITableHDU):
                hdu._container = self
        return has_new_hdu

    def get_version(self):
        return self._OI_VER

    def _verify(self, option='warn'):

        errors = super()._verify(option) 
       
        # check OI extensions are valid names and fit the OIFITS version
        for hdu in self[1:]:
            extname = hdu.header.get('EXTNAME', '')
            extrevn = hdu.header.get('OI_REVN', 0)
            if extname[0:3] == 'OI_' and not isinstance(hdu, _OITableHDU):
                err_text = f"Invalid OIFITS extention: {extname}"
                fix_text = "Replaced underscore by dash"
                def fix(hdu=hdu):
                    hdu.header['EXTNAME'] = 'OI-' + extname[3:]
                err = self.run_option(option, err_text=err_text,
                                  fix_text=fix_text, fix=fix)
                errors.append(err)
            if hdu._OI_VER != self._OI_VER: 
                err_text = f"Extension {extname} rev. {extrevn} in {clsname}"
                err = self.run_option(option, err_text=err_text, fixable=False) 
                errors.append(err)
       
        # silently fix EXTVER 
        self.update_extver()
        self.update_primary_header()

        return errors

    def get_HDUs(self, exttype, filter=None):

        hdus = [h for h in self[1:] if isinstance(h, exttype)]
        
        if filter is not None:
            hdus = [h for h in hdus if filter(h)]

        return hdus
    
    def get_dataHDUs(self):
        return self.get_HDUs(_DataHDU)

    def get_HDU(self, extype, filter=None):
        
        hdus = self.get_HDUs(extype, filter=filter)
        
        if not hdus:
            return None
        return hdus[0]

    def get_OITableHDUs(self):
        return self.get_HDUs(_OITableHDU)

    def get_arrayHDU(self, arrname):
        def same_arrname(h): return h.get_arrname() == arrname
        return self.get_HDU(_ArrayHDU, same_arrname)

    def get_targetHDU(self):
        return self.get_HDU(_TargetHDU) 

    def get_wavelengthHDUs(self):
        return self.get_HDUs(_WavelengthHDU)
    
    def get_wavelengthHDU(self, insname):
        def same_insname(h): return h.get_insname() == insname
        return self.get_HDU(_WavelengthHDU, same_insname)

    def get_corrHDU(self, corrname):
        def same_corrname(h): h.get_corrname() == corrname
        return self.get_HDU(_CorrHDU, same_corrname)

    def get_inspolHDU(self, arrname):
        def same_arrname(h): return h.get_arrname() == arrname
        return self.get_HDU(_InspolHDU, same_arrname)

    def to_table(self):
       
        tabs = [h._to_table(full_uv=True) for h in self.get_dataHDUs()]
        colnames = tabs[0].colnames 
        cols = [_ma.hstack([t[n] for t in tabs]) for n in colnames]

        tab = _table.Table(cols, names=colnames)

        for x in ['INT_TIME', 'U1COORD', 'U2COORD', 'V1COORD', 'V2COORD']:
            tab.columns[x].format = '7.3f'
        for x in ['EFF_WAVE', 'EFF_BAND']:
            tab.columns[x].format = '7.5e'
        tab.columns['MJD'].format = '7.5f' 
        for x in ['value', 'error']:
            tab.columns[x].format = '7.5g'

        return tab

    def append(self, h2):
        raise NotImplementedError()

    def copy(self):
        return  type(self)([h.copy() for h in self])

    def __add__(self, other):
        return merge(self, other)

    def get_referrers(self, hdu):

        refkey = getattr(hdu, '_REFERENCE_KEY', None)
        if refkey is None:
            return []
        
        refval = hdu.header[refkey]
        hdus = [h for h in self if h.header.get(refkey, '') == refval and h 
               is not hdu]
        return hdus
        

    def update_extver(self):

        extnames = _np.unique(h.header.get('EXTNAME', None) for h in self[1:])
        for extname in extnames:
            hdus = [h for h in self[1:] 
                            if h.header.get('EXTNAME', None) == extname]
            if len(hdus) > 1:
                for i, h in enumerate(hdus):
                    h.header['EXTVER'] = i + 1

    def update_primary_header(self):
        
        header = self[0].header
        datahdus = self.get_dataHDUs()
        
        # fix dates
        mjdobs = min(h.MJD.min() for h in datahdus)
        header['DATE-OBS'] = _Time(mjdobs, format='mjd').isot 
        
        # deal with keywords for atomic observations
        targets = self.get_targetHDU().data
        if len(targets) == 1:
            wavehdus = self.get_wavelengthHDUs()
            mjdend = max(h.MJD.max() for h in datahdus)
            header['OBJECT'] = targets['TARGET'][0] 
            header['RA'] = targets['RAEP0'][0]
            header['DEC'] = targets['DECEP0'][0]
            header['EQUINOX'] = targets['EQUINOX'][0]
            header['MJD-OBS'] = mjdobs
            header['MJD-END'] = mjdend
            header['WAVELMIN'] = min(h.EFF_WAVE.min() for h in wavehdus)
            header['WAVELMAX'] = max(h.EFF_WAVE.max() for h in wavehdus)
        else:
            for keyw in ['RA', 'DEC', 'UTC', 'LST', 'EQUINOX', 'RADECSYS', 
                'TEXPTIME', 'MJD-OBS', 'MJD-END', 'BASE_MIN', 'BASE_MAX', 
                'WAVELMIN', 'WAVELMAX', 'NUM_CHAN', 'VIS2ERR', 'VISPHERR', 
                'T3PHIERR']:
                if keyw in header:
                    del header[keyw]

        # Deal with all MULTI keywords
        if len(targets) > 1:
            header['OBJECT'] = 'MULTI'

        arrnames = _np.unique([h.get_arrname() for h in datahdus])
        if len(arrnames) == 1:
            header['TELESCOP'] = arrnames[0]
        else:
            header['TELESCOP'] = 'MULTI'

        insnames = _np.unique([h.get_insname() for h in datahdus])
        if len(insnames) == 1:
            header['INSMODE'] = insnames[0]
        else:
            header['INSMODE'] = 'MULTI'

            # deduce instrument
        ins = [_re.sub('([A-Za-z]+).*', '\\1', i).upper() for i in insnames]
        ins = _np.unique(ins)
        if len(ins) == 1:
            ins = ins[0]
        else:
            ins = 'MULTI'
        self[0].header['INSTRUME'] = ins

        # missing ones
        for keyw in ['PROG_ID', 'REFERENC', 'PROCSOFT', 'OBSTECH',
                     'OBSERVER', 'TELESCOP']:
            if keyw not in header:
                header[keyw] = 'UNKNOWN'


class OIFITS1(_OIFITS):
    _OI_VER = 1

class OIFITS2(_OIFITS):
    _OI_VER = 2

