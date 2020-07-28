from astropy.io.fits.hdu import base as _fitsbase
from .. import utils as _u

_InheritCardDescription = _u.InheritConstantArray(
                            '_CARDS',
                            dtype=[
                                ('name', 'U32'), ('required', bool),
                                ('test', object), ('default', object),
                            ]
                          )

class _ValidHDU(
    _InheritCardDescription,
    _fitsbase._ValidHDU
):
    
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
