from astropy.io.fits.hdu import base as _fitsbase
from astropy.io import fits as _fits
from .. import utils as _u
import copy as _copy

# All Valid HDUs for OIFITS will inherit a _CARDS structured array
# describing the specific FITS keywords in the header
_InheritCardDescription = _u.InheritConstantArray(
    '_CARDS',
    dtype=[
        ('name', 'U32'),     # keyword name
        ('required', bool),  # whether it must be present
        ('test', object),    # func that tests the validity of keyword value
        ('default', object), # default value. None: no default possible
        ('comment', 'U47')   # associated comment
    ]
)

class _ValidHDU(
    _InheritCardDescription,
    _fitsbase._ValidHDU
):

    def copy(self):
        """

Create a duplicate of an OIFITS extension or primary header, with data
and header are copied.

        """
        data = self.data
        if data is not None:
            data = data.copy()
        cp = type(self)(data=data, header=self.header.copy())
        if hasattr(self, '_container'):
            cp._container = self._container
        return cp

    def get_container(self):
        return getattr(self, '_container', None)


    def verify(self, option='warn'):

        with warnings.catch_warnings():
            warnings.simplefilter('always')
            super().verify(option=option)

    def _verify(self, option='warn'):

        errors = super()._verify(option)
       
        header = self.header
        for card in self._CARDS:
            name = card['name']
            if card['required'] or name in header:
                test = card['test']
                default = card['default']
                self.req_cards(name, None, test, default, option, errors)
        
        return errors

    def __ge__(self, other):

        return not self < other

    def __gt__(self, other):

        return not self <= other

    def __lt__(self, other):
    
        return self is not other and self <= other

    def __le__(self, other):

        # A OIFITS primary HDU comes first, then a FITS primary HDU,
        # then an OIFITS extension, and last other FITS extensions.
        if (isinstance(self, _fits.PrimaryHDU) or
            not isinstance(other, _ValidHDU)):
            return True
        
        if isinstance(other, _fits.PrimaryHDU):
            return False

        # OIFITS extensions will be ordered by type first
        order = {'OI_TARGET': 1, 'OI_ARRAY': 2, 'OI_WAVELENGTH': 3,
                 'OI_VIS': 4, 'OI_VIS2': 5, 'OI_T3': 6, 'OI_FLUX': 7,
                 'OI_CORR': 8, 'OI_INSPOL': 9}

        sheader = self.header
        oheader = other.header
        sext = sheader.get('EXTNAME', '')
        oext = oheader.get('EXTNAME', '')
        sorder = order.get(sext, 99999)
        oorder = order.get(oext, 99999)

        if sorder != oorder:
            return sorder < oorder

        # ordering uses reference names, date of observation and, in the
        # last resort, memory address.
        
        arr1 = sheader.get('ARRNAME', '')
        arr2 = sheader.get('ARRNAME', '')
        ins1 = sheader.get('INSNAME', '')
        ins2 = oheader.get('INSNAME', '')
        corr1 = sheader.get('CORRNAME', '')
        corr2 = oheader.get('CORRNAME', '')
        mjd1 = min(getattr(self.data, 'MJD', [0]))
        mjd2 = min(getattr(other.data, 'MJD', [0]))
        id1 = id(self)
        id2 = id(other)

        return (arr1, ins1, corr1, mjd1, id1) <  (arr2, ins2, corr2, mjd2, id2)

    def to_version(self, version):
       
        if self._OI_VER == version:
            return self

        newclass = self.get_class(version=version) 
        newobj =  newclass(data=self.data, header=self.header)

        return newobj 

    @classmethod
    def get_class(cls, version=2):
        """Get the class corresponding to OIFITS version"""
        
        if hasattr(cls, '_OI_VER'):
            cls = cls.__base__
        
        subclasses = cls.__subclasses__()
        for newcls in subclasses:
            if getattr(newcls, '_OI_VER', None) == version:
                break
        else:
            error = f"Could not find a class maching {version=}"
            raise RuntimeError(error)

        return newcls

class _OIFITS1HDU(_ValidHDU):
    _OI_VER = 1

class _OIFITS2HDU(_ValidHDU):
    _OI_VER = 2

