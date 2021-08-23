from astropy.io import fits
import numpy as np

from . import hdu
from . import fitsutils

from .hdu.base import _ValidHDU
from .hdu.table import _OITableHDU
from .hdu.data import _DataHDU
from .hdu.target import _TargetHDU, new_target_hdu, new_target_hdu_from_simbad
from .hdu.array import _ArrayHDU, new_array_hdu
from .hdu.t3 import _T3HDU, new_t3_hdu
from .hdu.vis2 import _Vis2HDU, new_vis2_hdu
from .hdu.vis import _VisHDU, new_vis_hdu
from .hdu.flux import _FluxHDU, new_flux_hdu
from .hdu.wavelength import _WavelengthHDU, new_wavelength_hdu
from .hdu.referenced import _Referenced
from .hdu.corr import _CorrHDU
from .hdu.inspol import _InspolHDU
from .hdu.primary import _PrimaryHDU, new_primary_hdu

__all__ = ["OIFITS1", "OIFITS2", "open", "openlist", "merge", 
           "set_merge_settings", "new_target_hdu", "new_array_hdu",
            "new_wavelength_hdu", "new_vis_hdu", "new_vis2_hdu",
            "new_t3_hdu", "new_flux_hdu", "new_target_hdu_from_simbad",
            "new_primary_hdu"]

def open(filename, mode='readonly', lazy_load_hdus=True, **kwargs):
    """

Open an OIFITS file.

Arguments
---------

filename (str)
    File name
mode (str, optional, default: 'readonly')
    Read mode
lazy_load_hdus (bool, optional, default: True)
    Whether to only load HDUs if needed

    """
    # default if cannot be determined
    cls = 'OIFITS2'
    
    # we need to read the primary HDU to see which FITS version it is
    if mode != 'ostream':
        with fits.open(filename, lazy_load_hdus=True) as hdulist:
            if list.__len__(hdulist): # len(h) would load all HDUs
                content = hdulist[0].header.get('CONTENT', '')
                cls = OIFITS2 if content == 'OIFITS2' else OIFITS1

    hdus = cls.fromfile(filename, lazy_load_hdus=lazy_load_hdus, **kwargs)

    return hdus

def openlist(filenames):
    """

Open a list of OIFITS files and merge them.

Arguments
---------

filenames (str × N)
    File names

    """
    hdulists = [open(f, lazy_load_hdus=False) for f in filenames]
    hdulist = _merge(*hdulists, _inplace=True)

    return hdulist

def merge(*hdulists):
    """

Merge several OIFITS.

Arguments
---------

hdulists1, hdulist2, ...
    OIFITS objects 

    """
    return _merge(*hdulists, _inplace=False)     

def _merge(*hdulists, _inplace=False):


    # Use the latest OIFITS version used by the OIFITS to be merged
    # We order them by newest version first.

    oiver = [getattr(hdulist, '_OI_VER', 0) for hdulist in hdulists]
    order = np.argsort(oiver)[::-1]
    hdulists = [hdulists[o] for o in order]
    
    
    maxver = oiver[order[0]]
    cls = type(hdulists[order[0]])

    # Copy files if necessary
    if not _inplace:
        hdulists = [hdulist.copy() for hdulist in hdulists]
    
    # Trick to avoid Delayed columns
    # for hdulist in hdulists:
    #    for hdu in hdulist:
    #        nope = hdu.data

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

    # * if some extensions are of the wrong version, convert them
    for i, hdu in enumerate(hdus):
        if isinstance(hdu, _ValidHDU):
            hdus[i] = hdu._to_version(maxver)
    
    # * update the primary HDU header inferring some keywords from 
    #   contents
    merged = cls(hdus)
    merged.update_primary_header()
    
    # * final tweaks
    merged.verify('silentfix+ignore')
    merged.sort()
    merged.update_extver()

    return merged

def _process_primaryHDUs(hdus):

    from .fitsutils import merge_fits_headers
    
    # merge primary headers of OIFITs primary HDUs 
    headers = [hdu.header for hdu in hdus if isinstance(hdu, _PrimaryHDU)]
    header = fitsutils.merge_fits_headers(*headers)
    hdus[0].header = header

    # delete void OIFIT Primary headers and convert others to images
    for i in range(len(hdus)-1, 0, -1):
        hdu = hdus[i]
        if isinstance(hdu, fits.PrimaryHDU):
            if hdu.data is not None:
                hdus[i] = fits.ImageHDU(hdu)
            else:
                del hdus[i]
       
def _remove_equal_OITableHDUs(hdus, cls=_OITableHDU):
    
    for i, hdu in enumerate(hdus):

        if not isinstance(hdu, cls):
            continue

        for j in range(len(hdus)-1, i, -1):
            if hdus[j] == hdu:
                del hdus[j]     


def _rename_conflicting_OITableHDUs(hdus, cls=_Referenced):

    from .utils import SequentialName

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
        
        refname1 = SequentialName(hdu1.header[refkey])
        refhdus2 = [hdu2 for hdu2 in hdus2 if hdu2 & hdu1]
        refnames2 = [SequentialName(h.header[refkey]) for h in refhdus2] 
        equal = [refname1 == refname2 for refname2 in refnames2]
        indices = np.argwhere(equal)[:,0]
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

            refhdu2 = refhdus2[index]
            new_refname2 = str(new_refnames2.pop(0))
            refhdu2.rename(new_refname2)

def _merge_OITableHDUs(hdus, cls=_OITableHDU):

    # find all equal HDUs and remove them

    for i, hdu in enumerate(hdus):    

        if not isinstance(hdu, cls):
            continue
        
        # Find all mergeable HDUs, merge them 
        nhdu = len(hdus)
        is_mergeable = [hdus[k] % hdu if k > i else False for k in range(nhdu)]
        where_mergeable = np.argwhere(is_mergeable)[:,0]
        if len(where_mergeable):
            mergeable = [hdus[m] for m in where_mergeable]
            hdus[i] = hdus[i].merge(*mergeable)
            for j in where_mergeable[::-1]:
                del hdus[j]
   
class _OIFITS(fits.HDUList):

    _merge_station_distance = 0.1    # 0.1 m
    _merge_array_distance = 10       # 10 m
    _merge_target_distance = 2.5e-8  # ~0.005 mas 
    _merge_target_name_match = False # allow different target designations.

    @staticmethod
    def set_merge_settings(*, station_distance=0.1, array_distance=10, 
           target_distance=2.5e-8, target_name_match=False):
        """

Set the default behaviour when merging several OIFITS objects

Arguments
---------

station_distance (float, default: 0.1)
    Maximum distance (m) for stations to be considered the same
array_distance (float, default: 10)
    Maxiumum distance (m) for array centres to be considered from
    the same array
target_distance (float, default: 2.5e-8 i.e. 5 mas)
    Maximum distance (deg) for two targets to be considered the
    same one
target_name_match (bool, default: False)
    Whether target names must match exactly to be considered the
    same one

        """
        _OIFITS._merge_station_distance = station_distance
        _OIFITS._merge_array_distance = array_distance
        _OIFITS._merge_target_distance = target_distance
        _OIFITS._merge_target_name_match = target_name_match

    def __repr__(self):

        name = type(self).__name__
        str_ = [str(h) for h  in self]
        return f"<{name} at {hex(id(self))}: {' '.join(str_)}>"

    def __str__(self):

        name = type(self).__name__
        str_ = [str(h) for h  in self]

        return f"<{name}: {' '.join(str_)}>"
   
    def __init__(self, hdus=[], file=None):
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
            if isinstance(hdu, (_OITableHDU, _PrimaryHDU)):
                hdu._container = self
        return has_new_hdu

    def get_version(self):
        return self._OI_VER

    def verify(self, option='warn'):

        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter('always')
            super().verify(option=option)

    def _verify(self, option='warn'):

        errors = super()._verify(option) 
      
        name = f"{type(self).__name__} ({hex(id(self))})"

        # check cross-references
        for hdu in [*self.get_dataHDUs(), *self.get_inspolHDUs()]:

            arrname = hdu.header.get('ARRNAME', None)
            if self.get_arrayHDU(arrname) is None:
                err_text = f'OI_ARRAY with ARRNAME={arrname} not found'
                err = sef.run_option(option, err_text=err_text, fixable=False)
                errors.append(err)

            insname = hdu.header.get('INSNAME', None)
            if self.get_wavelengthHDU(insname) is None:
                err_text = f'OI_WAVELENGTH with INSNAME={insname} not found'
                err = sef.run_option(option, err_text=err_text, fixable=False)
                errors.append(err)

        # check OI extensions are valid names and fit the OIFITS version
        for hdu in self[1:]:
            
            extname = hdu.header.get('EXTNAME', '')
            extrevn = hdu.header.get('OI_REVN', 0)
            
            if extname[0:3] == 'OI_' and not isinstance(hdu, _OITableHDU):
                err_text = f"Invalid OIFITS extention: {extname} in {name}"
                fix_text = "Replaced underscore by dash"
                def fix(hdu=hdu):
                    hdu.header['EXTNAME'] = 'OI-' + extname[3:]
                err = self.run_option(option, err_text=err_text,
                                  fix_text=fix_text, fix=fix)
                errors.append(err)
              
            if hdu._OI_VER != self._OI_VER: 
                err_text = f"Extension {extname} rev. {extrevn} in {name}"
                err = self.run_option(option, err_text=err_text, fixable=False) 
                errors.append(err)
       
        return errors

    def update_uv(self):
        """

Update the (u, v) coordinates (UCOORD, VCOORD) using information of
the array and target information contained in OI_ARRAY and OI_TARGET 
tables.

It is a low-precision routine not meant for high precision work.

Warnings
--------

(u, v) may differ from (UCOORD, VCOORD) determined by a data processing
software because

a. (UCOORD, VCOORD) can be averaged independently from MJD 
(see OIFITS standard)
b. Atmospheric refraction is dealt with approximately, while (UCOORD,
   VCOORD) may have none to full modelling of the atmosphere.

For the VLTI, the differences amount to 0.1-0.2% error on baselines
(a few centimetres).

        """
        for hdu in self.get_dataHDUs():
            hdu.update_uv()

    def get_HDUs(self, exttype, filter=None):
        """
Get all HDUs of a given extension type matching given criteria

Arguments
---------

exttype (type or str)
    Class of the extension or string with the extension name (e.g.
    OI_VIS2, OI_WAVELENGTH)
filter (func)
    Function taking an extension object and returning either True or 
    False

        """
        if isinstance(type(exttype), type):
            hdus = [h for h in self[1:] if isinstance(h, exttype)]
        elif isinstance(exttype, str):
            hdus = [h for h in self[1:] if h.header['extname'] == h]
        else:
            raise NotImplementedError('')

        if filter is not None:
            hdus = [h for h in hdus if filter(h)]

        return hdus
   
    def get_tableHDUs(self, filter=None):
        """
Get all HDUs containing an OI binary table
        """
        return self.get_HDUs(_OITableHDU, filter=filter)
        
 
    def get_dataHDUs(self, filter=None):
        """
Get all HDUs containing optical interferometry data.  

Keyword arguments
-----------------

filter
    Function taking an HDU and returning True/False.

        """
        return self.get_HDUs(_DataHDU, filter)

    def get_HDU(self, extype, /, filter=None):
        """
Get the first HDU of a given extension type matching given criteria

Arguments
---------

exttype (type)
    Extension type
filter (func)
    Function taking an extension object and returning either True or 
    False

        """
        hdus = self.get_HDUs(extype, filter=filter)
        
        if not hdus:
            return None
        return hdus[0]


    def get_arrayHDUs(self):
        return self.get_HDUs(_ArrayHDU)
    
    def get_arrayHDU(self, arrname):
        """
Get HDU containing array description (OI_ARRAY)

Arguments
---------

arrname (str)
    Name of the array

        """
        def same_arrname(h): return h.get_arrname() == arrname
        return self.get_HDU(_ArrayHDU, same_arrname)
    
    def get_vis2HDUs(self, filter=None):
        """
Get the HDUs containing square visibility amplitude data
        """
        return self.get_HDUs(_Vis2HDU, filter=filter)
    
    def get_visHDUs(self, filter=None):
        """
Get the HDUs containing square visibility data
        """
        return self.get_HDUs(_VisHDU, filter=filter)
    
    def get_t3HDUs(self, filter=None):
        """
Get the HDUs containing closure (3T) data
        """
        return self.get_HDUs(_T3HDU, filter=filter)

    def get_targetHDU(self):
        """
Get the HDU containing target information (OI_TARGET)
        """
        return self.get_HDU(_TargetHDU) 

    def get_fluxHDUs(self, filter=None):
        """
Get the HDUs containing flux information (OI_FLUX)
        """
        return self.get_HDUs(_FluxHDU)

    def get_wavelengthHDUs(self):
        """
Get all HDUs containing wavelength information (OI_WAVELENGTH)
        """
        return self.get_HDUs(_WavelengthHDU)

    def get_wavelengthHDU(self, insname):
        """
Get the HDU containing wavelength information (OI_WAVELENGTH) of a 
given instrumental setup.

Arguments
---------

insname (str)
    Name of the instrumental setup

        """ 
        def same_insname(h): return h.get_insname() == insname
        return self.get_HDU(_WavelengthHDU, same_insname)

    def get_corrHDUs(self):
        """
Get all HDUs containing correlation information (OI_CORR)
        """
        return self.get_HDUs(_CorrHDU)

    def get_corrHDU(self, corrname):
        """
Get the HDU containing a correlation matrix (OI_CORR)

Arguments
---------

corrname:
    Name of the correlation matrix

        """
        def same_corrname(h): h.get_corrname() == corrname
        return self.get_HDU(_CorrHDU, same_corrname)

    def get_inspolHDUs(self):

        return self.get_HDUs(_InspolHDU)

    def get_inspolHDU(self, arrname):
        """
Get the HDU containing the instrumental polarisation (OI_INSPOL)
corresponding to a given array.

Arguments
---------

arrname:
    Name of the array

        """
        def same_arrname(h): return h.get_arrname() == arrname
        return self.get_HDU(_InspolHDU, same_arrname)

    def _to_table(self, *, correlations=None, remove_masked=False,
        **kwargs):

        from astropy import table
        from scipy import sparse
        from numpy import ma 

        dataHDUs = self.get_dataHDUs()

        return_corr = correlations is not None

        # join all tables generated by each data HDU
        tabs = [h._to_table(full_uv=True, correlations=return_corr,
                    remove_masked=remove_masked, **kwargs) for h in dataHDUs]
        colnames = tabs[0].colnames 
        cols = [ma.hstack([t[n] for t in tabs]) for n in colnames]
        tab = table.Table(cols, names=colnames)
        for name in colnames:
            tab.columns[name].format = tabs[0].columns[name].format
        
        if not return_corr:
            return tab

        corr = sparse.identity(len(tab), format='dok') 

        # treat each OI_CORR separately, then look up indices
        # in the full table
        for corrHDU in np.unique(self.get_corrHDUs()):

            corrname = corrHDU.get_corrname()
            keep = tab['CORRNAME'] == corrname
            
            # global index: goes from 0 to len(tab) - 1
            # local index: CORRINDX tabulated for each OI_CORR
            gi = np.argwhere(keep)[:,0]
            li = tab['CORRINDX'][keep]

            xmatch = {l: g for l, g in zip(li, gi) if not li.mask}

            # looks quite pedestrian, can't I vectorise that?
            for line in corrHDU.data:
                i = xmatch[line['IINDX']]
                j = xmatch[line['JINDX']]
                val = line['CORR']
                if val != 0:
                    corr[i,j] = val
                    corr[j,i] = val

        tab.remove_columns(['CORRNAME', 'CORRINDX'])

        try:
            corr = getattr(corr, 'to' + correlations)()
        except AttributeError:
            raise ValueError(f"wrong correlation matrix format: {correlations}")

        return tab, corr

    def append(self, h2):
        """Not implemented"""
        raise NotImplementedError()

    def pop(self, index=-1):
        """Not implemented"""
        raise NotImplementedError()

    def copy(self):
        """

Create a duplicate of an OIFITS, without an attached file, with data 
and header are copied

        """
        return type(self)([h.copy() for h in self])

    def __add__(self, other):
        return merge(self, other)

    def trim(self, *, keep_ns_columns=False, 
            targets=None, target_filter=lambda targ: True,
            insnames=None, insname_filter=lambda ins: True,
            arrnames=None, arrname_filter=lambda arr: True,
            wavemin=0, wavemax=1e10, wave_filter=lambda wave: True):
        """

Trim unwanted wavelengths, instruments, instrumental setups, arrays, targets,
and non-standard columns from an OIFITS file

Arguments:
----------

keep_ns_columns (bool, default: False)
    Whether non-standard columns should be kept.  In the case they are
    kept, no trimming of wavelengths occurs if these columns have a spectral
    dimension.

wavemin (float, default: 0):
    Minimum wavelength to be kept

wavemax (float, default: 1e10):
    Maximum wavelength to be kept

wave_filter (func):
    A boolean function indicating whether a wavelength should be kept.  It
    may be combined with wavemin & wavemax.

targets (list of str):
    A list of targets to keep

target_filter (func):
    A boolean function indicating whether a target should be kept. It may be
    combined with targets.

arrnames (list of str):
    A list of array names to keep

arrname_filter (func):
    A boolean function indicating whether an array should be kept. It may be
    combined with arrnames
    
insnames (list of str):
    A list of instrument configurations to keep. 

insname_filter (func):
    A boolean function indicating whether an instrumental configuration 
    should be kept.  It may be combined with insnames.

Returns:
--------

    A trimmed OIFITS1 or OIFITS2

        """
        if any(isinstance(h, (_CorrHDU)) for h in self):
           msg = 'Trimming OIFITS with correlation info'
           raise NotImplementedError(msg)
        
        # Merge filter with other types of constraints

        def wfilter(wave):
            return (wave <= wavemax) & (wave >= wavemin) & wave_filter(wave)
        
        def ifilter(ins):
            keep = insname_filter(ins)
            if insnames is not None:
                keep &= ins in insnames
            return keep
        
        def afilter(arr):
            keep = arrname_filter(arr)
            if arrnames is not None:
                keep &= arr in arrnames
            return keep 
       
        def tfilter(targ):
            keep = target_filter(targ)
            if targets is not None:
                keep &= targ in targets
            return keep

        # Remove unwanted arrays

        trimmed = self # all HDUs are copied in later stage

        trimmed = [h for h in trimmed 
                if not (s := h.header.get('ARRNAME', '')) or afilter(s)]

        # Remove unwanted instrument configurations

        trimmed = [h for h in trimmed
                if not (s := h.header.get('INSNAME', '')) or ifilter(s)]

        # Remove unwanted wavelengths, targets, and insnames from 
        # lines and columns.   
 
        trimmed = [h._trim_helper(target_filter=tfilter, wave_filter=wfilter,
                insname_filter=ifilter, keep_ns_columns=keep_ns_columns)
            if isinstance(h, _OITableHDU) else h.copy() for h in trimmed]


        # Remove empty HDUs (either zero line or zero dimension in
        # wavelength-dependent data)

        trimmed = [h for h in trimmed 
            if not isinstance(h, _OITableHDU) or h.ndata()]
        
        return type(self)(trimmed)

    def bin_spectral_channels(self, R):
        """

Bin spectral channels down to a given spectral resolution.  HDUs with lower 
resolution are left untouched.   

If nchan is the number of resulting spectral channels, the resolution at the
middle of the band will be between R and R * (1 + 1 / nchan)

Argument
--------

R (float)
    Desired spectral resolution.  

Return value
------------

An OIFITS1 or OIFITS2 object.

        """
        if any(isinstance(h, (_CorrHDU, _InspolHDU)) for h in self):
           msg = 'Binning OIFITS with correlation or polarimetry info'
           raise NotImplementedError(msg)

        hdulist = [h.copy() for h in self 
                        if not isinstance(h, (_DataHDU, _WavelengthHDU))]

        for whdu in self.get_wavelengthHDUs():

            ins = whdu.get_insname()
            whdu, weights = whdu._bin_helper(R)
            
            dhdus = self.get_dataHDUs(filter=lambda h: h.get_insname() == ins)
            dhdus = [dhdu._bin_helper(weights) for dhdu in dhdus]

            hdulist += [whdu, *dhdus]

        return type(self)(hdulist)

    def update_extver(self):
        """

Update the EXTVER header keyword in extensions sharing the same name. 
While certainly unused by applications, this is a requirement of the
FITS standard.

        """
        extnames = np.unique(h.header.get('EXTNAME', None) for h in self[1:])
        for extname in extnames:
            hdus = [h for h in self[1:] 
                            if h.header.get('EXTNAME', None) == extname]
            if len(hdus) > 1:
                for i, h in enumerate(hdus):
                    h.header['EXTVER'] = i + 1

    def to_version(self, n):
        """
Convert an OIFITS object to version 1 or 2 of the standard.

Arguments
---------

n (int: 1 or 2)
    Version of the OIFITS standard

        """
        cls = type(self)
        for newcls in type(self).__base__.__subclasses__():
            if newcls._OI_VER == n:
                break
        else:
            raise TypeError(f'unknown OIFITS version: {n}')

        hdulist = newcls([hdu._to_version(n) for hdu in self])
        hdulist.verify('silentfix+ignore')
        hdulist.update_primary_header()

        return hdulist 

    def update_primary_header(self):
        """

Update the primary header to match the information in OI table 
extensions

        """ 
        self[0].update_header()

    def visualize(self, xvar, observable, /, *, 
            group_by=None, color_by=None, fig=None, **kwargs):
        """
Quick visualisation of an interferometric observable as a function of
time or baseline.  

Arguments
---------

xvar (str)
    x variable of the plot. It can be 'MJD', 'date', or 'time', 'baseline', 
    or 'spatial_frequency'.

observable (str)
    y variable. It can be any tabulated observable of the standard like
    VIS2DATA, T3AMP, FLUXDATA (OIFITS2), etc.

Keyword arguments
-----------------

fig (int or matplotlib.figure.Figure, optional)
    Matplotlib figure or figure number.  

group_by (str, optional)
    Subplots are grouped by 'target' or 'sta_config'.

color_by (str, optional)
    Data points are coloured according to 'target', 'sta_config', or
    'insname'.

target  (default: all are kept)
    Target name or list of target names

insname  (default: all are kept)
    Instrument configuration name or list thereof

arrname  (default: all are kept)
    Array name or list of array names

mjd_min  (default: all are kept)
    Minimum Modified Julian date

mjd_max  (default: all are kept)
    Maximum Modified Julian date

wavelmin (default: 0)
    Minimum wavelength

wavelmax (default: +inf)
    Maximum wavelength


Returns
-------

fig (matplotlib.figure.Figure)
    A matplotlib figure.

        """
        from matplotlib import pylab 

        res = self._plot_helper(xvar, observable, **kwargs)
        (target, inscfg, stacfg, mjd), x, (y, dy) = res

        if xvar in  ['MJD', 'mjd']:
            x = x.to_value('mjd')
        elif xvar in ['date', 'time']:
            x = x.to_datetime()

        if group_by is None:
            if xvar in ['MJD', 'mjd', 'date', 'time']:
                group_by = 'sta_config'
            else:
                group_by = 'target'

        if group_by == 'target':
            key = target
            subkeys = [stacfg, inscfg]
        else:
            key = stacfg
            subkeys = [target, inscfg]

        if color_by == 'target':
            subkey = target
        elif color_by == 'sta_config':
            subkey = stacfg
        elif color_by == 'insname':
            subkey = inscfg
        elif color_by == 'mjd':
            subkey = mjd
        else:
            if len(np.unique(subkeys[0])) > 1:
                subkey = subkeys[0]
            else:
                subkey = subkeys[1]

        unique_key, unique_index = np.unique(key, return_inverse=True)
        unique_subkey, unique_subindex = np.unique(subkey, return_inverse=True)
        
        if fig is None or isinstance(fig, int):
            fig = pylab.figure(fig)
            fig.clf()
        
        naxes = len(unique_key)
        ny = int(np.sqrt(2 * naxes))
        nx = int(np.ceil(naxes / ny))
        axes = fig.subplots(ny, nx, sharex=True, sharey=True, squeeze=False)
        
        fig.subplots_adjust(hspace=0, wspace=0)
        for j in range(0, ny):
            for i in range(0, nx):
                if j < nx - 1:
                    axes[j][i].set_xticklabels([])
                else:
                    axes[j][i].set_xlabel(xvar)
                if i > 0:
                    axes[j][i].set_yticklabels([]) 
                else:
                    axes[j][i].set_ylabel(observable)

        for i in range(len(unique_key)):
            k = unique_key[i]
            keep = unique_index == i
            ax = fig.axes[i]
            xi, yi, dyi = x[keep], y[keep], dy[keep]
            subindexi = unique_subindex[keep]
            for j in range(len(unique_subkey)):
                subk = unique_subkey[j]
                keep = subindexi == j
                xj, yj, dyj = xi[keep], yi[keep], dyi[keep]
                ax.errorbar(xj, yj, yerr=dyj, label=subk, linestyle='none')
            ax.text(0.50, 0.99, k, 
                ha='center', va='top', transform=ax.transAxes)
        ax.legend()

        return fig
 
    # this function returns the data to plot:
    # * x variable (time, projected baseline length, or spatial frequency)
    # * y variable (value and error of the observable)
    # * target
    # * baseline used in the format 'ARRAY/STATION1-...'
    def _plot_helper(self, xvar, observable,  /, **kwargs):

        from astropy.time import Time

        tab = self.to_table(remove_masked=True, observable=observable, **kwargs)

        if not tab:
            raise RuntimeError('no data match criteria')

        if xvar in ['MJD', 'date', 'time']:
            x = Time(tab['MJD'], format='mjd')
        elif xvar in ['baseline', 'spatial_frequency']:
            x = np.zeros_like(tab['MJD'])
            t3 = ~tab['U2COORD'].mask
            t2 = ~tab['U1COORD'].mask & ~t3 
            u1, v1 = tab['U1COORD'][t2], tab['V1COORD'][t2]
            b1 = np.sqrt(u1 ** 2 + v1 ** 2)
            x[t2] = b1
            u1, v1 = tab['U1COORD'][t3], tab['V1COORD'][t3]
            u2, v2 = tab['U2COORD'][t3], tab['V2COORD'][t3]
            u3, v3 = -u1-u2, -v1-v2
            b1 = np.sqrt(u1 ** 2 + v1 ** 2)
            b2 = np.sqrt(u2 ** 2 + v2 ** 2)
            b3 = np.sqrt(u3 ** 2 + v3 ** 2)
            x[t3] = abs(b1 * b2 * b3) ** (1/3)
            if xvar == 'spatial_frequency':
                x /= tab['EFF_WAVE']
        else:
            choices = 'MJD, date, time, baseline, or spatial_frequency'
            raise NotImplementedError(f"x variable must be {choices}")

        if not isinstance(x, Time):
            x = np.asarray(x)
        y = np.asarray(tab['value'])
        dy = np.asarray(tab['error'])
        
        arr = tab['ARRNAME'].tolist()
        cfg = tab['STA_CONFIG'].tolist()        
        stacfg = np.asarray([f"{a}/{b}" for a, b in zip(arr, cfg)])
        inscfg = np.asarray(tab['INSNAME'])        
        mjd = np.asarray(tab['MJD'])

        target = np.asarray(tab['TARGET'])

        return (target, inscfg, stacfg, mjd), x, (y, dy) 
 

class OIFITS1(_OIFITS):
    """Top-level class of Optical Interferometry FITS format, version 1."""
    _OI_VER = 1

    def insert_arrayHDU(self, ahdu, dataHDUs=None):
        """

Insert an OI_ARRAY extension to an OIFITS. 

Arguments:
----------

ahdu:
    OI_ARRAY extension

dataHDUs (list of HDUs):
    HDUs that will be referring to ahdu.  They must not previously be
    referring to an OI_ARRAY extension.

        """
        self.append(ahdu)
        arrname = ahdu.get_arrname()

        if dataHDUs is None:
            dataHDUs = self.get_dataHDUs()

        for dhdu in dataHDUs:

            if not isinstance(dhdu, _DataHDU):
                msg = 'dataHDUs must refer to OI_VIS, OI_VIS2, OI_T3 extensions'
            elif dhdu not in self:
                msg = 'dataHDUs must refer to extensions contained by OIFITS'
            elif dhdu.header.get('ARRNAME', None):
                msg = 'dataHDUs must not be already referring to an OI_ARRAY'
            else:
                dhdu.header['ARRNAME'] = arrname            
                continue
     
            raise ValueError(msg)

    def to_table(self, /, *, remove_masked=False, **kwargs):
        """

Convert to a flat table containing one scalar interferometric 
observable per line.  It is possible to only keep data corresponding
to a given list of targets, observables, instrumental configurations,
array configurations, wavelength range and/or date range. 


Arguments
---------

remove_masked (bool, default: False)
    Remove masked values.

observable (default: all are kept) 
    Observable or list of observables (e.g. VIS2DATA)

observable_type  (default: all are kept)
    Observable type or list of observable types (absolute, differential, etc.)

target  (default: all are kept)
    Target name or list of target names

insname  (default: all are kept)
    Instrument configuration name or list thereof

arrname  (default: all are kept)
    Array name or list of array names

mjd_min  (default: all are kept)
    Minimum Modified Julian date

mjd_max  (default: all are kept)
    Maximum Modified Julian date

wavelmin (default: 0)
    Minimum wavelength

wavelmax (default: +inf)
    Maximum wavelength


Returns
-------

tab (astropy.table.Table)
    A table with one scalar observable per line.

        """
        return self._to_table(remove_masked=remove_masked, **kwargs)


class OIFITS2(_OIFITS):
    """Top-level class of Optical Interferometry FITS format, version 2."""
    
    _OI_VER = 2

    def _verify(self, option='warn'):

        errors = super()._verify(option)

        # check cross-references
        for hdu in self.get_dataHDUs():

            corrname = hdu.header.get('CORRNAME', None)
            if corrname and self.get_corrHDU(corrname) is None:
                err_text = f'OI_CORR with CORRNAME={corrname} not found'
                err = sef.run_option(option, err_text=err_text, fixable=False)
                errors.append(err)

        return errors
    
    def to_table(self, /, *, correlations=None, remove_masked=False, **kwargs):
        """

Convert to a flat table containing one scalar interferometric 
observable per line. It is possible to only keep data corresponding
to a given list of targets, observables, instrumental configurations,
array configurations, wavelength range, and/or date range. 

Arguments
---------

correlations (bool, default: None)
    Type of correlation matrix
    * None: no matrix is returned
    * 'array': numpy array
    * 'dense': numpy dense matrix
    * 'dok': scipy.sparse matrix in Dictionary Of Keys format
    * 'csc': scipy.sparse matrix in Compressed Sparse Column format
    * 'csr': scipy.sparse matrix in Compressed Sparse Row format
    * 'coo': scipy.sparse matrix in COOrdinate format

remove_masked (bool, default: False)
    Remove masked values.

observable (default: all are kept) 
    Observable or list of observables (e.g. VIS2DATA)

observable_type  (default: all are kept)
    Observable type or list of observable types (absolute, differential, etc.)

target  (default: all are kept)
    Target name or list of target names

insname  (default: all are kept)
    Instrument configuration name or list thereof

arrname  (default: all are kept)
    Array name or list of array names

mjd_min  (default: all are kept)
    Minimum Modified Julian date

mjd_max  (default: all are kept)
    Maximum Modified Julian date

wavelmin (default: 0)
    Minimum wavelength

wavelmax (default: +inf)
    Maximum wavelength

Returns
-------

tab (astropy.table.Table)
    table with one scalar observable per line

corr (a matrix, optional)
    sparse correlation matrix 

        """
        return self._to_table(correlations=correlations, 
                            remove_masked=remove_masked, **kwargs)

set_merge_settings = _OIFITS.set_merge_settings

