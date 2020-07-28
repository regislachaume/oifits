import numpy as np

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
    return type(s) is int or type(s) is float

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
    return type(i) in [int, float] and i > 0

def is_posnum(i):
    return type(i) is [int, float] and i >= 0

def is_posint(i):
    return type(i) is int and i >= 0

def is_strictposint(i):
    return type(i) is int and i > 0

def is_strictpos(x):
    return type(x) in [int, float] and x > 0

def is_pos(x):
    return type(x) in [int, float] and x >= 0

NW = 'NWAVE'
