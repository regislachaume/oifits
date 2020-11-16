import numpy as np
from astropy.io import fits
import re
from scipy.spatial.transform import Rotation

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

def is_frame1(s):
    return s == 'GEOCENTRIC'
def is_frame2(s):
    return s in ['GEOCENTRIC', 'SKY']
def is_fovtype2(s):
    return s in ['FVWM', 'RADIUS']

def is_veltyp(s):
    return s in ['LSR', 'BARYCENT', 'HELIOCEN', 'GEOCENTR', 'TOPOCENT']
def is_veldef(s):
    return s in ['OPTICAL', 'RADIO']

def is_zero(i):
    return i == 0

def is_oifits2(s):
    return s in ['OIFITS2']

def is_oifits1(s):
    return s in ['OIFITS', 'OIFITS1']

def is_nonempty(s):
    return type(s) is str and len(s)

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

def is_num(s):
    return isinstance(s, (int, float, np.number))

def is_strictposnum(i):
    return is_num(i) and i > 0

def is_posnum(i):
    return is_num(i) and i >= 0

def is_pos(x):
    return is_num(x) and x >= 0

def is_strictpos(x):
    return is_num(x) and x > 0


def is_int(i):
    return isinstance(i, (int, np.integer))

def is_strictposint(i):
    return is_int(i) and i > 0

def is_posint(i):
    return is_int(i) and i >= 0


def is_category(s):
    return s in ['SCI', 'CAL']

def is_calstat(s):
    return s in ['U', 'C']

def is_fovtype(s):
    return s in ['RADIUS', 'FWHM']

def is_amptyp2(s):
    return s in ['absolute', 'differential', 'correlated flux']

def is_phityp2(s):
    return s in ['absolute', 'differential']


NW = 'NWAVE'

def store_default(columns, name, *, default, shape=()):
    if name not in columns:
        columns[name] = default
    if shape:
        columns[name] = _np.full(shape, columns[name])

def rotation3d(array, axes, angles, degrees=False):
    rot = Rotation.from_euler(axes, angles, degrees=degrees)
    shape = np.shape(array)
    array = rot.apply(array.reshape((-1, 3))).reshape(shape)
    return array
