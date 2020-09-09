import numpy as np
from astropy.io import fits
import re

class SequentialName(object):
    def __init__(self, s, n=1):
        s = str(s)
        if m := re.match("^(.*)_([0-9]+)$", s):
            s = m.groups()[0]
            n = int(m.groups()[1])
        self._s = s
        self._n = n
    def __str__(self):
        if self._n == 1:
            return self._s
        return f"{self._s}_{self._n}"
    def __repr__(self):
        cls = type(self).__name__
        return f"{cls}('{str(self)}')"
    def __int__(self):
        return self._n
    def __eq__(self, other):
        return self._s == other._s and self._n == other._n
    def __mod__(self, other):
        """Test if base names are equal"""
        return self._s == other._s
    def __add__(self, n):
        return type(self)(self._s, self._n + n)
    @classmethod
    def next_available(cls, used, n=None):
        if n is None:
            return cls.next_available(used, n=1)[0]
        m = set(range(1, len(used) + n + 1)) - set([int(u) for u in used])
        m = sorted(m)[0:n]
        return [cls(used[0], i) for i in m] 

def merge_fits_headers(*headers, req_keys=[]):
    
    header = headers[0].copy()

    for h in headers[1:]:
        for card in h.cards:
            (name, value, comment) = card
            if name in ['', 'HISTORY', 'COMMENT']:
                header.append(card)
            elif name not in header:
                if 'HISTORY' in header:
                    header.insert('HISTORY', card)
                else:
                    header.append(card)
            elif value != header[name]:
                if type(header[name]) == str:
                    header[name] = 'MULTI'
                elif name not in req_keys:
                    header.remove(name)

    return header

def rename_duplicates(hdulists, key):
    maps = [{} for h in hdulists]
    names = [[h.header[key] for h in hdulist] for hdulist in hdulists]
    flat_names = [y for x in names for y in x]
    for name in np.unique(flat_names):
        if sum(name == flat_names) > 2:
            k = 1
            for i, (names_i, hdulist) in enumerate(zip(names, hdulists)):
                j = np.argwhere([n == name for n in names_i])
                if len(j):
                    j = j[0][0]
                    name_i = f"{name}_{k:02}"
                    hdulist[j].header[key] = name_i
                    maps[i][name] = name_i
    return maps           


def InheritConstantArray(varname, dtype=None):
    """Create a class that makes a (constant) class array inheritable, i.e.
    derived classes' array will inherit of all parent classes' elements plus
    the ones they are explicitly given.

    InheritX = InheritConstantArray('X', value=['X', 'x'])
    class A(InheritX):
        X = ['A']  # will contain ['A', 'X', 'x']

    class B(InheritX):
        X = ['B'] # will contain ['B', 'X', 'x']

    class C(A, B):
        X = ['C'] # will contain ['C', 'A', 'B', 'X', 'x']

    """
    
    basename = f"_InheritConstantArray{varname}"
    
    def init_subclass(cls):

        # it seems __class__ won't work here.  I don't understand magic.
        allbases = cls.mro()
        base = [b for b in allbases if b.__name__ == basename][0]
        super(base, cls).__init_subclass__()

        # base._X will store original cls.X (or base._X_ <- cls._X)
        uvarname = f'_{varname}' if varname[0] != '_' else f'{varname}_'
        if not hasattr(base, uvarname):
            setattr(base, uvarname,  {base: []})
        stored_values = getattr(base, uvarname)
        stored_values[cls] = cls.__dict__.get(varname, [])

        # Then we set cls.X value from base classes
        bases = [b for b in allbases if issubclass(b, base)]
        values = [i for b in bases[::-1] for i in stored_values.get(b, [])]
        if dtype is not None:
            values = np.array(values, dtype=dtype)
        setattr(cls, varname, values)
   
    # dynamic type creation 
    dct = {'__init_subclass__': init_subclass}
    base = type(basename, (object,), dct)
    
    return base

def _is_frame1(s):
    return s == 'GEOCENTRIC'

def _is_frame2(s):
    return s in ['GEOCENTRIC', 'SKY']
def _is_fovtype2(s):
    return s in ['FVWM', 'RADIUS']

def is_zero(i):
    return i == 0

def is_oifits2(s):
    return s in ['OIFITS2']

def is_oifits1(s):
    return s in ['OIFITS', 'OIFITS1']

def is_nonempty(s):
    return type(s) is str and len(s)

def is_num(s):
    return isinstance(s, (int, float))

def is_one(s):
    return s == 1

def is_two(s):
    return s == 2

def is_date(s):
    return is_nonemptystr(s)

def is_amptyp_rev2(s):
    return s in ['absolute', 'differential', 'correlated flux']

def is_phityp_rev2(s):
    return s in ['absolute', 'differential']

def is_strictposnum(i):
    return isinstance(i, (int, float)) and i > 0

def is_posnum(i):
    return isinstance(i, (int, float)) and i >= 0

def is_posint(i):
    return isinstance(i, int) and i >= 0

def is_strictposint(i):
    return isinstance(i, int) and i > 0

def is_strictpos(x):
    return isinstance(x, (int, float)) and x > 0

def is_pos(x):
    return isinstance(x, (int, float)) and x >= 0

def is_category(s):
    return s in ['SCI', 'CAL']

NW = 'NWAVE'


def store_default(columns, name, *, default, shape):
    if name not in columns:
        columns[name] = default
    if _np.shape(columns[name]):
        columns[name] = _np.full(shape, columns[name])


def _get_fits_col_dtype(name, col):
    col = asarray(col)
    return (name, col.dtype.str, col.shape[1:])

fits_types = {
     'i2': '2-bytes int', 'i4': '4-bytes int',
     'f4': '4-bytes float', 'f8': '8-bytes float',
     'c8': '8-bytes complex', 'c16': '16-bytes complex',
     'i1': 'boolean', 'b1': 'boolean'
}
fits_formats = {
    'i2': 'I', 'i4': 'J',
    'f4': 'E', 'f8': 'D',
    'c8': 'C', 'c16': 'M',
    'i1': 'L', 'b1': 'L',
}


def dtype_descr(t):
    t = t.strip("|<>")
    if t[0] in 'SU':
        len_ = int(t[1:])
        type_ = 'string'
    elif m := re.match('(int|float|complex)([0-9]+)', t):
        len_ = int(m.groups()[1]) // 8
        type_ = m.groups()[0]
    else:
        return fits_types.get(t, '')
    descr =  f"{len_}-byte{'s' if len_ > 1 else ''} {type_}"
    return descr

def dtype_to_fits(t, shape):
    t = t.strip("|<>")
    if t[0] in 'SU':
        len_ = int(t[1:])
        fmt = 'A'
    else:
        len_ = 1
        fmt = fits_formats[t]
    return f"{len_ * int(_np.prod(shape[1:]))}{fmt}"

def fits_column(*, name=None, array=None, unit=None, shape=None, dtype=None,
                    null=None, **kwarg):

    if array is not None:
        array = np.asarray(array, dtype=dtype)
    if shape is None:
        shape = array.shape

    fmt = dtype_to_fits(array.dtype, shape)

    dim = None
    if len(shape) >= 3 - (fmt[-1] == 'A'):
        dim  = shape[1:]
        if fmt[-1] == 'A':
            dim = (int(fmt[:-1]), *dim)
        dim = str(dim)
        dim = re.sub(' ', '', dim)

    col = fits.Column(array=array, format=fmt, name=name, unit=unit, dim=dim,
                null=null)
    
    return col

