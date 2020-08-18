from .table import _OITableHDU, _OITableHDU1, _OITableHDU2
from .table import _OIFITS1HDU, _OIFITS2HDU
from .array import _MayHaveArrayHDU,_MustHaveArrayHDU
from .target import _MustHaveTargetHDU
from .wavelength import _MustHaveWavelengthHDU
from .corr import _MayHaveCorrHDU
from .inspol import _MayHaveInspolHDU

from .. import utils as _u

from astropy.io import fits
from astropy import table
from numpy import ma
import re as _re

import scipy
import copy
import numpy as np

class _DataHDU(_OITableHDU):
    
    _CARDS = [
        ('DATE-OBS', True, _u.is_nonempty, None),
    ]
    _COLUMNS = [
        ('TARGET_ID', True, '>i2', (),         _u.is_strictpos, None,  None), 
        ('MJD',       True, '>f8', (),         None,            None,  "d"),
        ('INT_TIME',  True, '>f8', (),         None,            None,  "s"), 
        ('FLAG',      True, '|b1', ('NWAVE',), None,            False, None),
    ]
    
    def get_obs_type(self, name, shape='data', flatten=False):

        return self._resize_data('N/A', shape, flatten)

    def get_uv(self, shape='table', flatten=False):

        uv = [self._resize_data(self.data[x], shape_flatten)
                    for x in self._COORD_COLUMNS]
    
    def _update_targetid(self, index_map):

        hdu = np.copy(self)
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
            x = ma.masked_array(x, mask=not hasattr(self, name))
            return x

        mask = self.FLAG
        if hasattr(self, name):
            x = getattr(self, name)
        else:
            x = self._resize_data(default, shape)
            flag = flag | True
        x = ma.masked_array(x, mask=mask)
        if flatten:
            x = x.ravel()
 
        return x 

    def get_reference_channels(self, shape='data', flatten=False):

        visref = self._resize_data(0, 'data', flatten)
        return ma.masked_array(visref, mask=True)

    def _to_table(self, full_uv=False):

        names = self._table_colnames(full_uv=full_uv)
        cols = self._table_cols(full_uv=full_uv)
        return table.Table(cols, names=names)

    @classmethod
    def _get_uvcoord_names(cls, full_uv=False):
       
        if full_uv:
            return ['U1COORD', 'V1COORD', 'U2COORD', 'V2COORD']
         
        names = [c for c in cls._COLUMNS['name'] if _re.match('[UV].?COORD', c)]
        names = sorted(names, key=lambda x: x[1::-1])
        return names

    def to_table(self):

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
        
        for ref in ['ARRNAME', 'INSNAME', 'CORRNAME', 'EXTNAME']:
            sref = self.header.get(ref, '')
            for other in others:
                oref = other.header.get(ref, '')
                if sref != oref:
                    txt = f'Cannot merge tables with {ref} = {sref} & {oref}'
                    raise RuntimeError(txt)

        return self._merge_helper(*others)
    
# OIFITS1 Table 
class _DataHDU1(
        _DataHDU, 
        _MustHaveTargetHDU,
        _MayHaveArrayHDU,
        _MustHaveWavelengthHDU,
        _OIFITS1HDU
      ):
    _COLUMNS = [('TIME', True, '>f8', (), None, 0., "s")]
        

# OIFITS1 Table rev1
class _DataHDU11( 
         _DataHDU1,
         _OITableHDU1,
      ):
    pass

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
    _COLUMNS = [('TIME', True, '>f8', (), _u.is_zero, 0., "s")]

