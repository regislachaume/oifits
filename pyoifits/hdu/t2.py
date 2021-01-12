from .data import _DataHDU
from .. import utils as _u

class _T2HDU(_DataHDU):
    
    _COLUMNS = [
        ('STA_INDEX', True, '2I', (2,), _u.is_strictpos, None, None,
            'station indices in matching OI_ARRAY table'),
        ('UCOORD',    True, '1D', (),   None,            None, "m",
            'u coordinate'),
        ('VCOORD',    True, '1D', (),   None,            None, "m",
            'v coordinate'),
    ]
    
    def __getattr__(self, name):

        if name == 'U1COORD':
            return self.data['UCOORD']
        if name == 'V1COORD':
            return self.data['VCOORD']
        
        return super().__getattr__(name)

    @classmethod
    def from_data(cls, *, insname, mjd, fits_keywords={}, **columns):
       
        _u.store_default(columns, 'ucoord', default=0.)
        _u.store_default(columns, 'vcoord', default=0.)
 
        return super().from_data(insname=insname, mjd=mjd,
            fits_keywords=fits_keywords, **columns)

    def get_uvw(self, refraction=False):
        """

Get the (u, v, w) coordinates of the baseline

Arguments
---------

refraction (bool, optional, default: False)
    Whether refraction is taken into account

Returns
-------

(u, v, w)
    Projected baseline in metres

Precision
---------
    A few centimetres for hectometric baselines, mainly due to the uncertainty
    on the time-averaging.
   
        """
        uvw = self.get_stauvw(refraction=refraction)
        uvw = uvw[:,1,:] - uvw[:,0,:]
        u, v, w = uvw.T

        return u, v, w

    def update_uv(self):

        if self.get_arrayHDU() is None:
            return

        u, v, w = self.get_uvw()
        self.UCOORD = u
        self.VCOORD = v
