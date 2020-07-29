from .base import _ValidHDU
from .. import utils as _u

from astropy.io import fits as _fits
from astropy import table
from numpy import ma

import re
import scipy
import copy
import numpy as _np 


def _get_fits_col_dtype(name, col):
    col = asarray(col)
    return (name, col.dtype.str, col.shape[1:])

_fits_types = {
     '>i2': '2-bytes int', '>i4': '4-bytes int',
     '>f4': '4-bytes float', '>f8': '8-bytes float',
     '>c8': '8-bytes complex', '>c16': '16-bytes complex',
     '|i1': 'boolean', '|b1': 'boolean'
}


def _dtype_descr(t):
    if t[1] in 'SU':
        size = int(t[2:])
        return f"{size}-byte{'s' if size > 1 else ''} string"
    descr = _fits_types.get(t, '')
    return descr

_InheritColumnDescription = _u.InheritConstantArray(
                              '_COLUMNS',
                              dtype=[
                                  ('name', 'U32'),  ('required', bool),
                                  ('type', object), ('shape', object),
                                  ('test', object), ('default', object), 
                                  ('unit', 'U16'),
                              ]
                            )

class _OITableHDU(
         _ValidHDU,
         _fits.BinTableHDU,
         _InheritColumnDescription,
      ):

    @classmethod
    def __init_subclass__(cls):
        super().__init_subclass__()
        if getattr(cls, '_EXTNAME', None) and  getattr(cls, '_OI_REVN', None): 
            _fits.hdu.base._BaseHDU.register_hdu(cls)
 
    @classmethod
    def match_header(cls, header):
        extname = getattr(cls, '_EXTNAME', None)
        oi_revn = getattr(cls, '_OI_REVN', None)
        if extname is None or oi_revn is None:
            return NotImplementedError
        return (_fits.BinTableHDU.match_header(header) and
                header.get('EXTNAME', '') == extname and
                header.get('OI_REVN', '') == oi_revn)

    def append_lines(self, lines):
      
 
        merged_data = _np.hstack([self.data, lines]) 
        merged = type(self)(data=merged_data, header=self.header)

        return merged

    #def __init__(self, data=None, header=None, name=None, container=None,
    #         cols=None, names=None, **kwargs):
    #   
    #    # A bit more flexibility to create from columns and/or data
    #    dtype = None
    #    if cols is not None:
    #        if isinstance(cols[0], fits.column.Column):
    #            names, cols = zip(*[(c.name, c.array) for c in cols])
    #        if isinstance(names, str):
    #            names = names.split(',')
    #        dtype = [_get_fits_col_dtype(n, c) for n, c in zip(names, cols)]
    #        data = list(zip(*cols))
    #    
    #    if len(data):
    #        data = scipy.rec.array(data, dtype=dtype)
    #        data = fits.BinTableHDU(data=data, **kwargs).data
    #    
    #    super().__init__(data=data, header=header, name=name, **kwargs)
    #    self.container = container

    def _xmatch(self, name, refhdu, refname, concatenate=False):
        """Helper to find target or array properties from indices"""
        
        ref_values = getattr(refhdu, name)
        ref_indices = getattr(refhdu, refname)
        xmatch = dict(zip(ref_indices, ref_values))
        indices = getattr(self, refname)
        shape = indices.shape + _np.shape(ref_values[0])
        values = _np.reshape([xmatch[i] for i in indices.flatten()], shape)
      
        if concatenate:
            values = [_np.atleast_1d(a) for a in values]
            values = ['-'.join([str(x) for x in a]) for a in values]

        return values

    def data_shape(self):

        nrows = len(self.data)
        
        # is there some wavelength data?
        colnames = self._get_spec_colnames()
        if not colnames:
            return (nrows,)

        return self.data[colnames[0]].shape

    def get_nwaves(self):

       return 0 

    def rename_columns(self, **names):
        
        self.data.dtype.names = [names.get(n, n) for n in self.data.dtype.names]
        for old, new in names.items():
            self.columns[old].name = new
        self.update()
 
    def _verify(self, option='warn'):
        
        errors = super()._verify(option)

        # Check all required columns are present
        for c in self._get_oi_colnames(required=True):
            if c not in self.data.names:
                name = self.__class__.__name__
                err_text = f"Missing column '{c}' in {name} object"
                err = self.run_option(option, err_text, fixable=False)
                errors.append(err)
           
        # Check column types
        for coldesc in self._get_oi_columns(required=False):
            name, req, type_, shape, test, default, unit = coldesc
            if name not in self.columns.names:
                continue
            
            # check the type
            dtype = self.data[name].dtype
            if type_ != dtype.str:
                spec = _dtype_descr(type_)
                real = _dtype_descr(dtype.str)
                err_text = f"Column '{name}' type must be {spec} but is {real}."
                fix_text = "Ignored."
                def fix(): pass
                err = self.run_option(option, err_text, fix_text, fix)
                errors.append(err)
            
            # check dimensionality
            dshape = self.data.dtype[name].shape
            if _u.NW in shape:
                nwave = self.get_nwaves()
                shape = tuple(nwave if d == _u.NW else d for d in shape)
            if shape != dshape:
                err_text = f"Column '{name}' should have dimension {shape} but has {dshape}."
                err = self.run_option(option, err_text, fixable=False) 
                errors.append(err)
            
            # Check values
            if test is not None:
                values = self.data[name]
                invalid = _np.array([not test(v) for v in _np.nditer(values)])
                invalid = invalid.reshape(values.shape)
                if invalid.any():
                    fixable = default is not None
                    if fixable:
                        try:
                            default = _np.array(default, dtype=dtype)
                            def fix(): values[invalid] = default
                        except:
                            fixable = False
                    if not fixable:
                        fix = None
                    val1 = values[invalid][0]
                    err_text = f"Column '{name}' has incorrect values. First encountered  '{val1}'."
                    fix_text = f"Replaced by default value"
                    err = self.run_option(option, err_text, fix_text, fix, fixable)
            
        # Non standard columns should start have prefix_
        colnames = self.data.dtype.names
        oi_colnames = self._get_oi_colnames()
        subst = {name: f"NS_{name}" for name in colnames 
                    if name not in oi_colnames and '_' not in name}
        if subst:
            nonstd = ', '.join([f"'{n}'" for n in subst.keys()])
            err_text = f"Column name(s) should start with prefix_ :'{nonstd}'."
            fix_text = "NS_ has been prefixed to column name(s)"
            def fix(h=self): h.rename_columns(**subst)
            err = self.run_option(option, err_text, fix_text, fix)
            errors.append(err)

        return errors

    def __repr__(self):

        return f"<{type(self).__name__} at {hex(id(self))} ({self._diminfo()})>"
   
    def __str__(self):
        
        return f"<{type(self).__name__} ({self._diminfo()})>"
 
    # Quick access to OICOLUMNS with hdu.VI2DATA, etc.
    def __getattr__(self, s):
      
        colnames = self._get_oi_colnames() 
        oicolnames = [x for x in self.columns.names if x in colnames]
        if s in oicolnames:
            return self.data[s][...]
        
        clsname = type(self).__name__
        err = f"'{clsname}' object has no attribute '{s}'"
        raise AttributeError(err)

    def zero(self):
        newhdu = self.copy()
        for name in self._get_spec_colnames():
            newhdu.data[name][...] = 0
        return newhdu
    
    def copy(self):
        cls = self.__class__
        result = cls.__new__(cls)
        for k, v in sorted(self.__dict__.items()):
            if k not in ['_container', '_file']:
                v = copy.deepcopy(v)
            setattr(result, k, v)
        return result
    
    def _diminfo(self):
        
        ncols = len(self.columns)
        nrows, *nwave = self.data_shape()
        if len(nwave):
            return f"{ncols}C×{nrows}R×{nwave[0]}W" 
        return f"{ncols}C×{nrows}R" 

    def _resize_data(self, x, shape='none', flatten=False, copy=True):
  
        data_shape = self.data_shape()
        if len(data_shape) == 1 or x is None:
            return x
        
        if shape == 'none':
            target_shape = ()
        elif shape == 'table':
            target_shape = (len(self.data),)
        elif shape == 'data':
            target_shape = data_shape
        else:
            raise ValueError(f"shape incorrect: '{shape}'")

        x_shape = _np.shape(x)
        if target_shape:
    
            if not x_shape: # scalar
                x = _np.full(target_shape, x)
            elif x_shape[0] == target_shape[0]:
                if len(target_shape) == 2:
                    x = _np.full((target_shape[1:] + x_shape), x).swapaxes(0,1)
            else:
                msg = f"dimension mismatch: {x_shape} and {target_shape}"
                raise ValueError(msg)

            if flatten:
                target_len = _np.prod(target_shape)
                flat_shape = (target_len, *x.shape[len(target_shape):])
                x = x.reshape(flat_shape)

        elif x_shape:

            x = x.copy()

        return x


    @classmethod
    def _get_oi_columns(cls, required=False, condition=None):
        cols = cls._COLUMNS
        if required:
            cols = cols[cols['required']]
        if condition is not None:
            cols = [c for c in cols if condition(c)]
        return cols
    
    @classmethod
    def _get_oi_colnames(cls, required=False, condition=None):
        cols = cls._get_oi_columns(required, condition)
        names = [c[0] for c in cols]
        return names
    
    @classmethod
    def _get_spec_columns(cls, required=False):
        return cls._get_oi_columns(required, lambda c: c[3] == ('NWAVE',))
    
    @classmethod
    def get_error_names(cls, required=False):
        cols = cls._get_spec_colnames(required)
        return [c for c in cols if c[-3:] == 'ERR']

    @classmethod
    def get_observable_names(cls, required=False):
        cols = cls._get_spec_colnames(required)
        return [c for c in cols if c[-3:] != 'ERR' and c != 'FLAG']
    
    @classmethod
    def _get_spec_colnames(cls, required=False):
        cols = cls._get_spec_columns(required)
        return [c[0] for c in cols]

    def __add__(self, other):

        header = self.header.copy()
        data1 = self.data

        # receive lines or HDU? 
        if isinstance(hdu, _fits.BinTableHDU):
            data2 = other.data
            header2 = other.header
        else:
            data2 = other
            header2 = None
        
        # merge headers if relevant
        if header2 is not None:
            for k,v,c in zip(header2, header2.values(), header.comments):
                if k not in header:
                    header.append(k, v, c)
        
        # merge columns and create zero-filled FITS table
        nrows1 = len(data1)  
        nrows = nrows1 + len(data2)
        merged_colnames = np.unique(data1.names, data2.names)
        cols = [data1.columns[c] if c in data1.names else data2.columns[c] 
                                            for c in merge_colnames]
        merged = self.from_columns(cols, nrows=nrows, header=header, fill=True)

        # copy data
        for name in data1.names:
            merged.data[name][:nrows1] = data1[name]
        for name in data2.names:
            merged.data[name][nrows1:] = data2[name]
        
        return merged
    
    def _merge(self, other, id_key, eq_keys, dist_keys, max_distance=0):

        # Keep the IDs in the first merged table and find available new
        # IDs in the second one
        old_id1 = getattr(self, id_key)
        nrows1 = len(old_id1)
        old_id2 = getattr(other, id_key)
        candidate_id2 = set(range(len(old_id1, old_id2))) - set(old_id1)
        candidate_id2 = sorted(list(candidate_id2))

        # Find if any row of the second table matches the first
        #  * it must exactly match all eq_keys (e.g. TARGET)
        #  * it must be close to coordinates in dist_keys (e.g. RA, DEC)
        new_id2 = _np.zeros_like(old_id2)
        kept_lines = [] 
        i = 0
        for j, row2 in enumerate(other.data):
            eq = [self.data[k] == row2[k] for k in eq_keys]
            coo = [(self.data[k] - row2[k]).T for k in dist_keys]
            coo = np.reshape(coo, (-1, nrows1))
            same = _np.logical_and(eq, axis=0)
            close = _np.linalg.norm(coo, axis=0) <= max_distance
            where = np.argwhere(same*close)
            if len(where): # if a match no new row
                new_id2[j] = old_id1[where[0,0]]
            else:
                new_id2[j] = candidate_id2[i]
                i += 1

        # merge tables
        kept_lines = _np.array(kept_lines)
        merged = self._append_lines(other.data[kept_lines])
        merged.data[id_key][nrows1:] = new_id2[kept_lines]

        # maps old to new indices for the second table
        map2 = {o: n for o, n in zip(old_id2, new_id2)}

        return merged, map2



class _OITableHDU1(_OITableHDU):
    _OI_REVN = 1
    _CARDS = [('OI_REVN', True, _u.is_one, 1)]

class _OITableHDU2(_OITableHDU):
    _OI_REVN = 2
    _CARDS = [('OI_REVN', True, _u.is_two, 2)]

_InitialisedLater = None

class _OIFITS1HDU(_OITableHDU):
    _OI_VER = 1

class _OIFITS2HDU(_OITableHDU):
    _OI_VER = 2

class _OITableHDU11(
        _OIFITS1HDU,
        _OITableHDU1
      ):
    pass

class _OITableHDU21(
        _OIFITS2HDU,
        _OITableHDU1
      ):
    pass

class _OITableHDU22(
        _OIFITS2HDU,
        _OITableHDU2
      ):
    pass
