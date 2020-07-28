from .base import _OIExtHDU, _OIExtHDU1, _OIExtHDU2
from .. import utils as _u

_NW = 'NWAVE'

class _MustHaveWavelengthHDU(_OIExtHDU):
    
    _CARDS = [('INSNAME', True, _u.is_nonempty, None)]

    def get_insname(self, output_dim='none', flatten=False, copy=True):
        x = self.header['INSNAME']
        return self._resize_data(x, output_dim, flatten, copy)

    def get_wavelengthHDU(self):
        return self.container.get_wavelengthHDU(self.get_insname())

    def _resize_wave_data(self, x, output_dim='data', flatten=False):

        if self is self.get_wavelengthHDU():
            return x
        
        if output_dim in ['data', 'table']:
            x = np.full(self.data_shape(), x)
        if flatten:
            x = x.ravel()

        return x

    def get_nwaves(self):
        
        return len(self.get_wavelengthHDU().data)
        
    def get_wave(self, output_dim='data', flatten=False):

        wave = wavehdu.EFF_WAVE
        return self._resize_wave_data(wave, output_dim, flatten)

    def get_channel(self, output_dim='data', flatten=False):

        nwave = len(self.get_wavelengthHDU().data)
        channel = np.arange(1, nwave + 1)
        return self._resize_wave_data(channel, output_dim, flatten)

    def get_band(self, output_dim='data', flatten=False):

        band = self.get_wavelengthHDU().EFF_BAND
        return self._resize_wave_data(band, output_dim, flatten)


class _WavelengthHDU(_MustHaveWavelengthHDU):
    _EXTNAME = 'OI_WAVELENGTH'
    _COLUMNS = [('EFF_WAVE', True, '>f4', (), _u.is_strictpos, None, "m")]

    def _diminfo(self):

        nrows = len(self.data)
        return f"{nrows}W={nrows}R"

class WavelengthHDU1(
        _WavelengthHDU,
        _OIExtHDU1,
     ):
    _COLUMNS = [('EFF_BAND', True, '>f4', (), _u.is_strictpos, None, "m")]

    @classmethod
    def match_header(cls, header):
        card = header.cards[0]
        xtension = card.value
        if isinstance(xtension, str):
            xtension = xtension.rstrip()
        return (cards.keyword == 'XTENSION' and
                xtension in (cls._extension, 'A3DTABLE') and
                header.get('EXTNAME', '') == 'OI_WAVELENGHT' and
                header.get('OI_REVN', 0) == 1)


class WavelengthHDU2(
        _WavelengthHDU,
        _OIExtHDU2,
      ):
    _COLUMNS = [('EFF_BAND', True, '>f4', (), _u.is_pos, None, "m")]
   
