from astropy.io import fits
from astropy import table
from numpy import ma

import re
import scipy
import copy
import numpy as np

from .. import utils as _u

function = type(lambda x: 0)

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
    if t[1] == 'S':
        size = int(t[2:])
        return f"{size}-byte{'s' if size > 1 else ''} string"
    descr = _fits_types.get(t, '')
    return descr


def _xmatch(hdu, name, refhdu, refname, concatenate=False):

    values = getattr(refhdu, name)
    xmatch = dict(zip(getattr(refhdu, refname), values))
    ids = getattr(hdu, refname)
    shape = ids.shape + np.shape(values[0])
    val = np.reshape([xmatch[i] for i in ids.flatten()], shape)
  
    if concatenate:
        val = [np.atleast_1d(a) for a in val]
        val = ['-'.join([str(x) for x in a]) for a in val]

    return val

_InheritCardDescription = _u.InheritConstantArray(
                            '_CARDS', 
                            dtype=[
                                ('name', 'U32'), ('required', bool),
                                ('test', object), ('default', object),
                            ]
                          )

_InheritColumnDescription = _u.InheritConstantArray(
                              '_COLUMNS',
                              dtype=[
                                  ('name', 'U32'),  ('required', bool),
                                  ('type', object), ('shape', object),
                                  ('test', object), ('default', object), 
                                  ('unit', 'U16'),
                              ]
                            )

class _OITableHDU(fits.BinTableHDU):
    
    @classmethod
    def match_header(cls, header):
        if not super().match_header(header):
            return False
        if not hasattr(cls, '_EXTNAME') or not hasattr(cls, '_OI_REVN'):
            return False
        return (header.get('EXTNAME', '') == cls._EXTNAME and
                header.get('OI_REVN', 0) == cls._OI_REVN)
    
    def __init__(self, data=None, header=None, name=None, container=None,
             cols=None, names=None, **kwargs):
       
        # A bit more flexibility to create from columns and/or data
        dtype = None
        if cols is not None:
            if isinstance(cols[0], fits.column.Column):
                names, cols = zip(*[(c.name, c.array) for c in cols])
            if isinstance(names, str):
                names = names.split(',')
            dtype = [_get_fits_col_dtype(n, c) for n, c in zip(names, cols)]
            data = list(zip(*cols))
        
        if len(data):
            data = scipy.rec.array(data, dtype=dtype)
            data = fits.BinTableHDU(data=data, **kwargs).data
        
        super().__init__(data=data, header=header, name=name, **kwargs)
        self.container = container

    def data_shape(self):

        nrows = len(self.data)
        
        # is there some wavelength data?
        colnames = self._get_spec_colnames()
        if not colnames:
            return (nrows,)

        return self.data[colnames[0]].shape

    def get_nwaves(self):

       return 0 

    def _chk_card(self, option, name, text, errors, exist=True, remove=False):
    
        if (name in self.header) is exist:
            return 

        err_text = 'missing' if exist else 'unexpected'
        err_text = f"{text}: keyword '{name}' is {err_text}"
        fix_text = "card removed"
        fixable = not exist and remove and name in self.header
        def fix(h=self.header): h.remove(name)
        error =self.run_option(option, err_text, fix_text, fix, fixable)
        errors.append(error)

    def _chk_col(self, option, name, text, errors, exist=True):

        if (name in self.columns.names) is exist:
            return

        err_text = 'missing' if exist else 'unexpected'
        err_text = f"{text}: column '{name}' is {err_text}"
        error = self.run_option(option, err_text, fixable=False)
        errors.append(error)

    def rename_columns(self, **names):
        
        self.data.dtype.names = [names.get(n, n) for n in self.data.dtype.names]
        for old, new in names.items():
            self.columns[old].name = new
        self.update()
 
    def _verify(self, option='warn'):
        
        errors = super()._verify(option)

        # Check revision number
        revn_ok = lambda r: r == self._REVN
        self.req_cards('OI_REVN', None, revn_ok, self._REVN, option, errors)
        
        # Check all required columns are present
        for c in self._get_oi_colnames(required=True):
            if c not in self.data.names:
                name = self.__class__.__name__
                err_text = f"Missing column '{c}' in {name} object"
                err = self.run_option(option, err_text, fixable=False)
                errors.append(err)
           
        # Check all required cards are present are have correct type
        for card in self._CARDS:
            name = card['name']
            value = self.header.get(name, None)
            if card['required'] or value is not None:
                test, default = card['test', 'default']
                fix = None
                if default:
                    def fix(h=self.header): h[name] = default
                self.req_cards(name, None, test, fix, option, errors) 

        # Check column types
        dtypes = {c[0]: c[1] for c in self._get_oi_columns()}
        shapes = {c[0]: c[2] for c in self._get_oi_columns()}
        tests = {c[0]: c[3] for c in self._get_oi_columns()}
        for name, type, *dim in self.data.dtype.descr:
            # only test OI columns
            if name not in dtypes:
                continue 
            # Check the types
            if type != dtypes[name]:
                spec = _dtype_descr(dtypes[name])
                real = _dtype_descr(type)
                err_text = f"Column '{name}' type should be {spec} but is {real}."
                fix_text = "Ignored."
                def fix(): pass
                err = self.run_option(option, err_text, fix_text, fix)
                errors.append(err)
            # Has the data the expected dimension.  'NWAVE' string
            # must be replace by number of wavelengths in INSNAME HDU
            dim = dim[0] if len(dim) else ()
            nwave = self.get_nwaves()
            ref_dim = tuple(nwave if d == 'NWAVE' else d for d in shapes[name])
            if dim != ref_dim:
                err_text = f"Column '{name}' should have dimension {ref_dim}."
                err = self.run_option(option, err_text, fixable=False) 
            # Check values
            test = tests[name]
            if test is None:
                continue
            values = self.data[name].ravel()
            invalid = np.array([not test(v) for v in values])
            if any(invalid):
                val1 = values[invalid][0]
                err_text = f"Column '{name}' has incorrect values. First encountered  '{val1}'."
                err = self.run_option(option, err_text, fixable=False)
            
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

        name = re.sub('_rev[12]', '', type(self).__name__)
        addr = hex(id(self))
        info = self._diminfo()
        return f"<{name} ({info})>"
    
    # Quick access to OICOLUMNS with hdu.VI2DATA, etc.
    def __getattr__(self, s):
      
        colnames = self._get_oi_colnames() 
        oicolnames = [x for x in self.columns.names if x in colnames]
        if s in oicolnames:
            return self.data[s][...]
        
        cls = type(self).__name__
        err = f"'type(self).__name__' object has no attribute '{s}'"
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
            if k not in ['container', '_file']:
                v = copy.deepcopy(v)
            setattr(result, k, v)
        return result
    
    def _diminfo(self):
        
        shape = self.data_shape()
        if len(shape) == 1:
            return f"{shape[0]}R" 
        return f"{shape[1]}W×{shape[0]}R"

    def _resize_data(self, x, output_dim='none', flatten=False, copy=True):
  
        shape = self.data_shape()
        if len(shape) == 1 or x is None:
            return x
        
        if output_dim == 'none':
            output_dim = ()
        elif output_dim == 'table':
            output_dim = (len(self.data),)
        elif output_dim == 'data':
            output_dim = shape
        else:
            raise ValueError(f"output_dim incorrect: '{output_dim}'")

        x_dim = np.shape(x)
        if output_dim:
    
            if not x_dim: # scalar
                x = np.full(output_dim, x)
            elif x_dim[0] == output_dim[0]:
                if len(output_dim) == 2:
                    x = np.full((output_dim[1:] + x_dim), x).swapaxes(0, 1)
            else:
                msg = f"dimension mismatch: {x_dim} and {output_dim}"
                raise ValueError(msg)

            if flatten:
                target_len = np.prod(output_dim)
                flatdim = (target_len, *x.shape[len(output_dim):])
                x = x.reshape(flatdim)

        elif x_dim:

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
        return cls._get_oi_columns(required, lambda c: c[2] == ('NWAVE',))
    
    @classmethod
    def get_error_names(cls, required=False):
        cols = cls._get_spec_colnames(cls, required)
        return [c for c in cols if c[-3:] == 'ERR']

    @classmethod
    def get_observable_names(cls, required=False):
        cols = cls._get_spec_colnames(cls, required)
        return [c for c in cols if c[-3:] != 'ERR' and c != 'FLAG']
    
    @classmethod
    def _get_spec_colnames(cls, required=False):
        cols = cls._get_spec_columns(required)
        return [c[0] for c in cols]

    def __add__(self, hdu2):

        nrows1 = len(self.data)
        if isinstance(hdu, fits.BinTableHDU):
            data2 = hdu.data2
        else:
            data2 = hdu
        nrows = nrows1 + len(data2)
        merged = self.from_columns(self.columns, nrows=nrows)
        for colname in self.columns.names:
            merged[colname][nrows1:] = data2[colname]
        
        return merged

class _OIExtHDU1(_OIExtHDU):
    _OI_REVN = 1
    _CARDS = [('OI_REVN', True, _u.is_one, 1)]

class _OIExtHDU2(_OIExtHDU):
    _OI_REVN = 2
    _CARDS = [('OI_REVN', True, _u.is_two, 2)]
