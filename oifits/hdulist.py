from astropy.io import fits as _fits
from astropy import table as _table
from numpy import ma as _ma
import numpy as _np

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
from .hdu.primary import PrimaryHDU1, PrimaryHDU2

from .hdu.data import _DataHDU

def oifitsopen(filename, mode='readonly', **kwargs):
    if mode == 'ostream':
        cls = OIFITS2
    else:
        with _fits.open(filename, mode=mode, lazy_load_hdus=True) as h:
            if list.__len__(h): # len(h) loads all HDUs
                content = h[0].header.get('CONTENT', '')
                cls = OIFITS2 if content == 'OIFITS2' else OIFTS1
            else:
                cls = OIFITS2
    return cls.fromfile(filename, mode=mode, **kwargs)

class _OIFITS(_fits.HDUList):
    """Top-level class of Optical Interferometry FITS format.
When a fits file is opened a HDUList object is returned."""

    def __repr__(self):

        name = type(self).__name__
        str_ = [str(h) for h  in self]
        return f"<{name} at {hex(id(self))}: {' '.join(str_)}>"

    def __str__(self):

        name = type(self).__name__
        str_ = [str(h) for h  in self]
        return f"<{name}: {primary} {' '.join(str_)}>"
   
    #def __init__(self, hdus=[], file=None, *, _copy_hdus=True):
    #  
    #    new_hdus = [hdus[0]] 
    #    cls = type(self)
    #    for h in hdus[1:]:
    #        extname = h.header.get('EXTNAME', '')
    #        if extname in self._OI_EXT:
    #            oiext = self._OI_EXT[extname]
    #            if _copy_hdus:
    #                h = oiext(data=h.data, header=h.header, _container=self)
    #            else:
    #                h.__class__ = oiext
    #                h._container = self
    #        new_hdus.append(h)
    #    
    #    _fits.HDUList.__init__(self, hdus=new_hdus, file=file)
    #

    # original _read_next_hdu() uses super().append(), ruining any clean 
    # attempt to subclass HDUList
    def _read_next_hdu(self):
        
        has_new_hdu = super()._read_next_hdu()
        if has_new_hdu:
            last_index = list.__len__(self) - 1 # len(x) will load all HDUs
            hdu = self[last_index]
            if isinstance(hdu, _OITableHDU):
                hdu._container = self
        return has_new_hdu

    def append(self, hdu):

        if type(hdu) in self._OI_EXT.values():
            hdu._container = self
        super().append(hdu)

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

        return errors

    def get_HDUs(self, exttype, filter=None):

        if type(exttype) is str:
            exttype = self._OI_EXT[extype]
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
    
    def get_arrayHDU(self, arrname):
        def same_arrname(h): return h.get_arrname() == arrname
        return self.get_HDU(_ArrayHDU, same_arrname)

    def get_targetHDU(self):
        return self.get_HDU(_TargetHDU) 

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

    #def __add__(self, hdulist2):
    #
    #    hdulist1 = self
    #    
    #    # target ID
    #    target1 = hdulist1.get_targetHDU()
    #    target2 = hdulist2.get_targetHDU()
    #    target, map1, map2 = target1._merge(target2)
    # 
    #    # data HDUs
    #    data1 = [h._update_targetid(map1) for h in hdulist1.get_targetHDUs()]
    #    data2 = [h._update_targetid(map2) for h in hdulist2.get_targetHDUs()]
    #
    #    # wavelengths HDUs
    #    for wave2 in hdulist2.get_wavelengthHDUs():
    #        pass 

    def update_extver(self):

        extnames = _np.unique(h.header.get('EXTNAME', None) for h in self[1:])
        for extname in extnames:
            hdus = [h for h in self[1:] 
                            if h.header.get('EXTNAME', None) == extname]
            if len(hdus) > 1:
                for i, h in enumerate(hdus):
                    h.header['EXTVER'] = i + 1

class OIFITS1(_OIFITS):
    _OI_VER = 1

class OIFITS2(_OIFITS):
    _OI_VER = 2

