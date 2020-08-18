from .table import _OITableHDU, _OITableHDU11, _OITableHDU22
from .. import utils as _u
import numpy as _np

_NW = 'NWAVE'

class _MustHaveWavelengthHDU(_OITableHDU):
    
    _CARDS = [('INSNAME', True, _u.is_nonempty, None)]

    def get_insname(self, shape='none', flatten=False, copy=True):
        x = self.header['INSNAME']
        return self._resize_data(x, shape, flatten, copy)

    def get_wavelengthHDU(self):
        return self._container.get_wavelengthHDU(self.get_insname())

    def _resize_wave_data(self, x, shape='data', flatten=False):

        if self is self.get_wavelengthHDU():
            return x
        
        if shape in ['data', 'table']:
            x = _np.full(self.data_shape(), x)
        if flatten:
            x = x.ravel()

        return x

    def get_nwaves(self):
        
        return len(self.get_wavelengthHDU().data)
        
    def get_wave(self, shape='data', flatten=False):

        wave = self.get_wavelengthHDU().EFF_WAVE
        return self._resize_wave_data(wave, shape, flatten)

    def get_channel(self, shape='data', flatten=False):

        nwave = len(self.get_wavelengthHDU().data)
        channel = _np.arange(1, nwave + 1)
        return self._resize_wave_data(channel, shape, flatten)

    def get_band(self, shape='data', flatten=False):

        band = self.get_wavelengthHDU().EFF_BAND
        return self._resize_wave_data(band, shape, flatten)


class _WavelengthHDU(_MustHaveWavelengthHDU):
    
    _EXTNAME = 'OI_WAVELENGTH'
    _REFERENCE_KEY = 'INSNAME'
    _COLUMNS = [('EFF_WAVE', True, '>f4', (), _u.is_strictpos, None, "m")]
    
    def _diminfo(self):

        nrows = len(self.data)
        return f"{super()._diminfo()}={nrows}W"

    def is_referred_to_by(self, other):
        return (not isinstance(other, _WavelengthHDU) and
                isinstance(other, _MustHaveWavelengthHDU) and
                self.get_insname() == other.get_insname())

    def __mod__(self, other):

        return False

class WavelengthHDU1(
        _WavelengthHDU,
        _OITableHDU11,
     ):
    _COLUMNS = [('EFF_BAND', True, '>f4', (), _u.is_strictpos, None, "m")]


class WavelengthHDU2(
        _WavelengthHDU,
        _OITableHDU22,
      ):
    _COLUMNS = [('EFF_BAND', True, '>f4', (), _u.is_pos, None, "m")]
   
