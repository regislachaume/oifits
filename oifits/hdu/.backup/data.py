from .base import _OIExtHDU, _OIExtHDU1, _OIExtHDU2
from .array import _MayHaveArrayHDU,_MustHaveArrayHDU
from .target import _MustHaveTargetHDU
from .wavelength import _MustHaveWavelengthHDU
from .corr import _MayHaveCorrHDU
from .inspol import _MayHaveInspolHDU

from .. import utils as _u

from astropy.io import fits
from astropy import table
from numpy import ma

import scipy
import copy
import numpy as np

class _DataHDU(_OIExtHDU):
    
    _CARDS = [
        ('DATE-OBS', True, _u.is_nonempty, None),
    ]
    _COLUMNS = [
        ('TARGET_ID', True, '>i2', (),         _u.is_strictpos, None,  None), 
        ('MJD',       True, '>f8', (),         None,            None,  "d"),
        ('INT_TIME',  True, '>f8', (),         None,            None,  "s"), 
        ('FLAG',      True, '|i1', ('NWAVE',), None,            False, None),
    ]
    
    def get_obs_type(self, name, output_dim='data', flatten=False):

        return self._resize_data('N/A', output_dim, flatten)

    def get_uv(self, output_dim='table', flatten=False):

        uv = [self._resize_data(self.data[x], output_dim_flatten)
                    for x in self._COORD_COLUMNS]
    
    def _update_targetid(self, index_map):

        hdu = np.copy(self)
        hdu.TARGET_ID = [index_map[id] for id in hdu.TARGET_ID]
        return hdu
    
    def _table_colnames(self, full_uv=False): 

        COORD = self._COORD_COLUMNS
        if full_uv:
            COORD = ['U1COORD', 'V1COORD', 'U2COORD', 'V2COORD']
            
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
        def hstack(x): return ma.hstack(x)
        
        for name in names:
            if name == 'value':
                col = hstack([getf(n) for n in obs_names])
            elif name == 'error':
                col = hstack([getf(n) for n in err_names])
            elif name == 'observable':
                col = np.hstack([resize(n) for n in obs_names])
            elif name == 'type':
                col = np.hstack([gett(n) for n in obs_names])
            else:
                col = hstack([getf(name)] * len(obs_names))
            cols.append(col)

        return cols

    def get_field(self, name, output_dim='none', flatten=False, default=None):

        if name == 'INSNAME':
            return self.get_insname(output_dim, flatten)
        if name == 'ARRNAME':
            return self.get_arrname(output_dim, flatten, default)
        if name == 'CORRNAME':
            return self.get_corrname(output_dim, flatten, default)
        if name == 'TARGET':
            return self.get_target(output_dim, flatten, default)
        if name == 'EFF_WAVE':
            return self.get_wave(output_dim, flatten)
        if name == 'EFF_BAND':
            return self.get_band(output_dim, flatten)
        if name == 'CHANNEL':
            return self.get_channel(output_dim, flatten)
        if name == 'REF_CHANNEL_BITFIELD':
            return self.get_reference_channels(output_dim, flatten)
        if name == 'STA_CONFIG':
            return self.get_sta_config(output_dim, flatten)

        DATACOLS = self._get_spec_colnames()
        DATACOLS.remove('FLAG')

        if name not in DATACOLS:
            if hasattr(self, name):
                x = getattr(self, name)
            else:
                x = default
            x = self._resize_data(x, output_dim, flatten)
            x = ma.masked_array(x, mask=not hasattr(self, name))
            return x

        mask = self.FLAG
        if hasattr(self, name):
            x = getattr(self, name)
        else:
            x = self._resize_data(default, output_dim)
            flag = flag | True
        x = ma.masked_array(x, mask=mask)
        if flatten:
            x = x.ravel()
 
        return x 

    def get_reference_channels(self, output_dim='data', flatten=False):

        visref = self._resize_data(0, 'data', flatten)
        return ma.masked_array(visref, mask=True)

    def _to_table(self, full_uv=False):

        names = self._table_colnames(full_uv=full_uv)
        cols = self._table_cols(full_uv=full_uv)
        return table.Table(cols, names=names)

    def to_table(self):

        tab = self._to_table()

        for x in ['INT_TIME', 'U1COORD', 'U2COORD', 'V1COORD', 'V2COORD',
                    'UCOORD', 'VCOORD']:
            if x in tab.colnames:
                tab.columns[x].format = '7.3f'
        for x in ['EFF_WAVE', 'EFF_BAND']:
            tab.columns[x].format = '7.5e'
        tab.columns['MJD'].format = '7.5f'
        for x in ['value', 'error']:
            tab.columns[x].format = '7.5g'
        
        return tab
    
    def _verify(self, option='warn'):
        
        errors = super()._verify(option)

        # STA_INDEX must be >= 0 
        sta_index = getattr(self, 'STA_INDEX', None)
        if sta_index is not None:
            sta_index = np.unique(sta_index)
            if any(sta_index <= 0):
                err_text = "'STA_INDEX' must be ≥ 0"
                err = self.run_option(option, err_text, fixable=False)     
            
        # STA_INDEX must be referenced
        refhdu = self.get_arrayHDU()
        if refhdu not in [None, self] and sta_index is not None:
            for i in np.unique(sta_index):
                if i not in refhdu.STA_INDEX:
                    err_text = "'STA_INDEX' not referenced in ArrayHDU: {i}"
                    self.run_option(option, err_text, fixable=False)
                    errors.append(err_text)
        
        return errors

class _DataHDU1(
        _DataHDU, 
        _MustHaveTargetHDU,
        _MayHaveArrayHDU,
        _MustHaveWavelengthHDU,
      ):
    _COLUMNS = [('TIME', True, '>f4', (), None, 0., "s")]
        

# OIFITS1 Table rev1
class _DataHDU11( 
         _DataHDU1,
         _OIExtHDU1,
      ):
    pass

class _DataHDU2(
        _DataHDU,
        _MustHaveTargetHDU,
        _MustHaveArrayHDU,
        _MustHaveWavelengthHDU,
        _MayHaveCorrHDU,
        _MayHaveInspolHDU 
      ):
    pass

# OIFITS2 Table rev1 (new table in OIFITS2)
class _DataHDU21(
        _DataHDU2,
        _OIExtHDU1
      ):
    pass

# OIFITS2 Table rev2 (table was available in OIFITS1 but was updated)
class _DataHDU22(
        _DataHDU2,
        _OIExtHDU2
      ):
    _COLUMNS = [('TIME', True, '>f4', (), _u.is_zero, 0., "s")]
