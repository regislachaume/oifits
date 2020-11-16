from .. import utils as _u
from .. import fitsutils as _fu

from .base import _ValidHDU, _OIFITS1HDU, _OIFITS2HDU

from astropy.io import fits as _fits
import numpy as _np 
import re as _re

# All OIFITS tables will inherit a _COLUMNS structured array describing
# the columns specified in the standard
_InheritColumnDescription = _u.InheritConstantArray(
    '_COLUMNS',
    dtype=[
        ('name', 'U32'),     # TTYPEn (column name, e.g. VIS2DATA)
        ('required', bool),  # whether the column must be present
        ('format', 'U3'),    # TFORMATn (data format, e.g. 1J or 15A)
        ('shape', object),   # a tuple including 'NWAVE' if necessary
        ('test', object),    # a func that tests the validity of the value
        ('default', object), # default value, if None: no default possible
        ('unit', object),    # TUNITn. Notes:
                             # * None: no unit or empty value for TUNITn
                             # * 'any': TUNITn may be present with any value
                             # * a value: TUNITn must be present with this value
        ('comment', 'U47'),  # Comment for TTYPEn in the header
    ]
)


class _OITableHDU(
         _ValidHDU,
         _fits.BinTableHDU,
         _InheritColumnDescription,
      ):

    _REFERENCE_KEY = None

    def __init__(self, data=None, header=None, uint=False, ver=None, 
                    character_as_bytes=False):
        _fits.BinTableHDU.__init__(self, data=data, header=header, uint=uint,
                    ver=ver, character_as_bytes=character_as_bytes)

    def update(self):
        
        super().update()

        # update the comments for the standard OIFITS columns and other
        # keywords
        header = self.header
        header.set('EXTNAME', self._EXTNAME, 'OIFITS extension name')

        columns = self._get_oi_columns()
        comments = header.comments

        for index, name in enumerate(self.columns.names, start=1):
           
            ttype = f"TTYPE{index}"
            if not comments[ttype]: 
                comment = f"name of column {index}"
                if name in columns['name']:
                    column = columns[columns['name'] == name][0]
                    comment = column['comment']
                comments[ttype] = comment
           
            tform = f"TFORM{index}"
            if not comments[tform]:
                comments[tform] = f"format for {name}"
            
            tunit = f"TUNIT{index}"
            if tunit in header and not comments[tunit]: 
                comments[tunit] = f"unit of {name}"
            
            tdim = f"TDIM{index}"
            if tdim in header and not comments[tdim]:
                comments[tdim] = f"dimension of {name}"
       
        # standard v.2 makes this optional, but it's good practice
        # isn't it? 
        self.add_datasum()
        self.add_checksum()

    @classmethod
    def from_columns(cls, columns, header=None, nrows=0, fill=False,
            character_as_bytes=False):

        columns = cls._fix_column_types(columns)

        return super().from_columns(columns, header=header, nrows=nrows,
            fill=False, character_as_bytes=character_as_bytes) 

    def to_version(self, version):
        """

Transform an OIFITS table between versions of the OIFITS standard.

Arguments
---------

version (int in 1..2)
    Version

Returns
-------

Binary table of same EXTNAME conforming to OIFITS standard #version.

Warning
-------

Losses may happen because
1) some column widths are shorter in v. 1 resulting in truncation.
2) some OIFITS2 columns do not exist in OIFITS1 and will be 
    prefixed with NS_ if transforming 2 -> 1.  The 1 -> 2 
    transform will not restore the original column names.

        """

        newobj = super().to_version(version)
        if hasattr(self, '_container'):
            newobj._container = self._container
        newobj.fix_column_types()

        return newobj

    @classmethod
    def __init_subclass__(cls):

        super().__init_subclass__()
        if getattr(cls, '_EXTNAME', None) and  getattr(cls, '_OI_REVN', None): 
            _fits.hdu.base._BaseHDU.register_hdu(cls)
 
    @classmethod
    def match_header(cls, header):
        """

Internal astropy.io.fits cuisine that was inadvertently exposed by the
project.

        """
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
    
    def __and__(self, other):

        return getattr(other, '_EXTNAME', '') == getattr(self, '_EXTNAME', None)

    def __mod__(self, other):

        h1, h2 = self.header, other.header

        return (self & other and
                h1.get('INSNAME', '') == h2.get('INSNAME', '') and
                h1.get('ARRNAME', '') == h2.get('ARRNAME', '') and
                h1.get('CORRNAME', '') == h2.get('CORRNAME', ''))
 
    def _xmatch(self, refhdu, refname, *, name=None, concatenate=False):
        """Helper to find target or array properties from indices"""
       
        ref_indices = getattr(refhdu, refname)
        if name is None:
            ref_values = list(range(len(ref_indices))) 
        else:
            ref_values = getattr(refhdu, name)
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
        """

Return the number of spectral channels

        """ 
        return 0 

    def rename_columns(self, **names):
        """

Rename several columns at once.

Syntax: tab.rename_columns(oldname1=newname1, ...)

        """
        self.data.dtype.names = [names.get(n, n) for n in self.data.dtype.names]
        oi_colnames = self._get_oi_columns(required=True)['name']

        for old, new in names.items():

            if old in oi_colnames:
                err_txt = 'Cannot rename a mandatory OIFITS column'
                raise RuntimeError(err_txt)
            if new in self.columns.names:
                err_txt = 'Cannot rename column to an existing one'
                raise RuntimeError(err_txt)

            self.columns[old].name = new

        self.update()
 
    def _verify(self, option='warn'):
        
        errors = super()._verify(option)

        OI_COLUMNS = self._get_oi_columns(required=False)

        # Check all required columns are present
        for c in self._get_oi_colnames(required=True):
            if c not in self.data.names:
                name = self.__class__.__name__
                err_text = f"Missing column '{c}' in {name} object"
                err = self.run_option(option, err_text, fixable=False)
                errors.append(err)
           
        # Check column types
        for coldesc in OI_COLUMNS:
            name, req, fmt, shape, test, default, unit, comment = coldesc
            if name not in self.columns.names:
                continue
            
            # check the type
            real = self.columns[name].format
            if len(fmt) == 1:
                real = real[-1]
            if fmt != real:
                err_txt = f"Column {name}: type must be {fmt} but is {real}."
                fix_txt = "Will try to fix."
                def fix(): pass
                err = self.run_option(option, err_txt, fix_txt, fix)
                errors.append(err)

            # check unit
            # * 'any': unit can be absent or present with any value
            # * None: no unit, so TUNITn may be absent or empty string
            # * other value: must be matched in TUNITn
            real = self.columns[name].unit
            spec = unit
            if (spec is None and real not in [None, ''] or
                spec not in [None, 'any'] and spec != real): 
                err_txt = f"Column {name}: unit must be {spec} but is {real}."
                fix_txt = "Fixed." 
                def fix(col=self.columns[name]): col.unit = spec
                err = self.run_option(option, err_txt, fix_txt, fix)
                errors.append(err)
            
            # check dimensionality
            dshape = self.data.dtype[name].shape
            if _u.NW in shape:
                nwave = self.get_nwaves()
                shape = tuple(nwave if d == _u.NW else d for d in shape)
            if shape != dshape:
                err_txt = f"Column {name}: shape is {dshape}, should be {shape}"
                err = self.run_option(option, err_txt, fixable=False) 
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
                    err_txt = f"Column '{name}' has incorrect values. "
                    err_txt += f"First encountered  '{val1}'."
                    fix_txt = f"Replaced by default value"
                    err = self.run_option(option, err_txt, fix_txt, fix, fixable)

        # Non standard columns should start have prefix_
        colnames = self.data.dtype.names
        oi_colnames = self._get_oi_colnames()
        subst = {name: f"NS_{name}" for name in colnames 
                    if name not in oi_colnames and name[0:3] != 'NS_'}
        if subst:
            nonstd = ', '.join([f"'{n}'" for n in subst.keys()])
            err_text = f"Column name(s) should start with prefix_ :'{nonstd}'."
            fix_text = "NS_ has been prefixed to column name(s)"
            def fix(h=self): h.rename_columns(**subst)
            err = self.run_option(option, err_text, fix_text, fix)
            errors.append(err)

        # Fix column types (string length, float/double)
        if 'fix' in option:
            self.fix_column_types()

        return errors

    @classmethod
    def _fix_column_types(cls, columns):

        oi_columns = cls._get_oi_columns(required=False)
        new_columns = []

        for col in columns:

            coldesc = oi_columns[oi_columns['name'] == col.name]

            if len(coldesc):
                coldesc = coldesc[0]
                try:
                    col = _fu.ascolumn(col, name=coldesc['name'],
                        unit=coldesc['unit'], format=coldesc['format'])
                except:
                    pass
            new_columns.append(col)

        return new_columns

    def fix_column_types(self):
        """

Fix the column data types if float or string widths do not conform with
the standard.

        """
        columns = self._fix_column_types(self.columns)
        fixed = any(a is not b for a, b in zip(columns, self.columns))
        
        if fixed: 
            # astropy.io.fits is a mess.
            #  both .data and .columns must be accessed.
            self.data = _fits.FITS_rec.from_columns(columns)
            for col, format in zip(self.columns, self.data.formats):
                col.format = format
            self.update()

    def __repr__(self):

        return f"<{type(self).__name__} at {hex(id(self))} ({self._diminfo()})>"
   
    def __str__(self):
        
        return f"<{type(self).__name__} ({self._diminfo()})>"
 
    # Quick access to OICOLUMNS with hdu.VI2DATA, etc.
    def __getattr__(self, s):
      
        oicolnames = self._get_oi_colnames() 
        if s in oicolnames and s in self.columns.names:
            return self.data[s][...]
        
        clsname = type(self).__name__
        err = f"'{clsname}' object has no attribute '{s}'"
        raise AttributeError(err)

    def __setattr__(self, s, v):

        oicolnames = self._get_oi_colnames()
        if s in oicolnames and s in self.columns.names:
            self.data[s] = v
        else:
            self.__dict__[s] = v

    def zero(self):

        newhdu = self.copy()
        for name in self._get_spec_colnames():
            newhdu.data[name][...] = 0
        return newhdu
    
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

    def merge(self, *others):
        return self._merge_helper(*others)

    def _merge_helper(self, *others, id_name=None, equality=lambda a, b: None):
        """Merge a set of OIFITS tables of the same kind.  id_name: column
ID that must be kept unique. equality: criteria to discard redundant rows.
        """

        # Check we are merging the same type of extension 
        ext1 =  self.header['EXTNAME']
        for other in others:
            ext2 = other.header['EXTNAME']
            if ext1 != ext2:
                txt = 'cannot merge FITS extensions {ext1} with {ext2}' 
                raise TypeError(txt)

        # Determine class when merging different revisions of a table
        hdus = [self, *others]
        i_maxrevn = _np.argmax([x._OI_REVN for x in hdus])
        cls = type(hdus[i_maxrevn])

        # Merged tables will keep all columns.  Values will be zero if not 
        # defined in one of the tables
        colnames, columns = _fu.merge_columns(*hdus)

        # Merge sets of FITS rows.
        # * In each set, rows duplicating one of the previous set is eliminated
        # * For each set a map of old_id -> new_id is built to avoid
        #   duplicate IDs.
        rows = [hdu.data for hdu in hdus]
        rows, maps = _fu.merge_rows(*rows, id_name=id_name, equality=equality)
        nrows = sum(len(r) for r in rows)

        # Merge headers.  
        headers = [hdu.header for hdu in hdus]
        req_keys = cls._CARDS['name'][cls._CARDS['required']]
        header = _fu.merge_fits_headers(*headers, req_keys=req_keys)
        
        # Create an empty merged fits with the right number of rows,
        # then fill it.
        #
        # FIXME null values!
        #
        
        merged = cls.from_columns(columns, nrows=nrows, fill=1, header=header)
        
        rowmin = 0
        for i, (data, map) in enumerate(zip(rows, maps)):
        
            ndata = len(data)
            rowmax = rowmin + ndata

            for name in data.names:
                if map and name == id_name:
                    values = [map.get(i, i) for i in data[name]]
                else:
                    values = data[name]
                merged.data[name][rowmin:rowmax][...] = values 
            rowmin += ndata

        # Update in HDUs refering to other
        for hdu, map in zip(hdus, maps):
            if map and (container := getattr(hdu, '_container', None)):
                for h in container.get_OITableHDUs():
                    if h.refers_to(hdu): 
                        field = h.data[id_name]
                        for old, new in map.items():
                            field[field == old] = new
        
        return merged 

    def refers_to(self, other):
        """Whether an OIFITS table makes a reference to another one via
        ARRNAME, INSNAME, or CORRNAME"""
        return other.is_referred_to_by(self)
    
    def is_referred_to_by(self, other):
        """Whether an OIFITS table is a reference for another table via
        ARRNAME, INSNAME, or CORRNAME"""
        return False

    def __eq__(self, other):

        if not isinstance(other, _OITableHDU):
            return False

        header1, header2 = self.header, other.header
        data1, data2 = self.data, other.data

        return (header1['EXTNAME'] == header2['EXTNAME'] and
                header1.get('INSNAME', '') == header2.get('INSNAME', '') and
                header1.get('ARRNAME', '') == header2.get('ARRNAME', '') and
                header1.get('CORRNAME', '') == header2.get('CORRNAME', '') and
                len(data1) == len(data2) and 
                data1.dtype == data2.dtype and
                (data1 == data2).all())

    @classmethod
    def _get_column_shape(cls, **columns):
        
        for col in cls.get_oi_columns(required=True):
            if ((name := col['name'] in columns) and
                (shape := col['shape']) and
                (value := columns[name]) is not None):
                shape = _np.shape(value)
                return shape

        raise RuntimeError('Cannot determine table shape from input')

    @classmethod
    def _guess_shape(cls, columns):

        obs_names = cls.get_observable_names()
        if len(obs_names):
            for obs_name in obs_names:
                if obs_name in columns:
                    shape = _np.shape(columns[obs_name])
                    return shape
            raise RuntimeError('cannot guess shape from data')
        
        oi_columns = cls._get_oi_columns()
        for name in oi_columns['name']:
            if name in columns:
                shape = _np.shape(columns[name])
                return shape
        raise RuntimeError('cannot guess shape from data')
         
    @classmethod
    def from_data(cls, *, version=None, fits_keywords={}, **columns):

        if not hasattr(cls, '_OI_VER'):
            if version is None:
                version = getattr(cls, '_OI_VER', 2)
            cls = cls.get_class(version=version)

        # FITS keywords and column names are upper case 
        fits_keywords = {k.upper(): v for k, v in fits_keywords.items()
                                                    if v is not None}
        columns = {k.upper(): v for k, v in columns.items()}
       
        # prefix non-standard columns
        oi_colnames = cls._get_oi_columns(required=False)['name']
        columns = {n if n in oi_colnames else f"NS_{n}": v
                    for n, v in columns.items()}
 
        # Header
        header = _fits.Header()
        for card in cls._CARDS:
            name = card['name']
            comment = card['comment']
            if name in fits_keywords:
                header.set(name, fits_keywords[name], comment)
                del fits_keywords[name]
            elif card['required']:
                if card['default'] is not None:
                    header.set(name, card['default'], comment)
        for name, value in fits_keywords.items():
            if not isinstance(value, (tuple, list)):
                value = [value]
            header.set(name, *value)
    
        # Guess shape
        shape = cls._guess_shape(columns)
        nrows = shape[0]
        
        def full(s, x):
            if _np.ndim(x) > len(s):
                return _np.asarray(x)
            return _np.full(s, x)

        # Guess errors from data    
        obs_names = cls.get_observable_names() 
        err_names = cls.get_error_names()
        for obs, err in zip(obs_names, err_names):
            if obs in columns:
                if err not in columns:
                    columns[err] = 0. 
        fcols = []

            # official OIFITS columns
        for col in cls._get_oi_columns(required=False):
            name = col['name']
            if name not in columns:
                if col['required']:
                    columns[name] = col['default']
                else:
                    continue
            if name in [*obs_names, *err_names, 'FLAG']:
                new_shape = shape
            else:
                new_shape = (nrows,)
            array = full(new_shape, columns[name]) 
            del columns[name]
            fcol = _fu.ascolumn(array, format=col['format'], unit=col['unit'],
                name=name)
            fcols.append(fcol)

            # additional columns
        for name, array in columns.items():
            fcol = _fu.ascolumn(array, name=name)
            fcols.append(fcol)
        
        tab = super().from_columns(fcols, header=header)
        return tab
    

class _OITableHDU1(_OITableHDU):
    _OI_REVN = 1
    _CARDS = [('OI_REVN', True, _u.is_one, 1, 
        '1st revision of this table in OIFITS format')]
    
class _OITableHDU2(_OITableHDU):
    _OI_REVN = 2
    _CARDS = [('OI_REVN', True, _u.is_two, 2, 
        '2nd revision of this table in OIFITS format')]

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
