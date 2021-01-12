"""
Implementation of the OI_WAVELENGTH binary table extension containing
the instrumental spectral setup.
"""

from .table import _OITableHDU, _OITableHDU11, _OITableHDU22
from .referenced import _Referenced

from .. import utils as _u
import numpy as _np
from astropy.io import fits as _fits


_NW = 'NWAVE'

class _MustHaveWavelengthHDU(_OITableHDU):
    
    _CARDS = [('INSNAME', True, _u.is_nonempty, None, 
        'instrumental setup name for cross-reference')]

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


class _WavelengthHDU(_MustHaveWavelengthHDU,_Referenced):
    
    _EXTNAME = 'OI_WAVELENGTH'
    _REFERENCE_KEY = 'INSNAME'
    _COLUMNS = [
        ('EFF_WAVE', True, '1E', (), _u.is_strictpos, None, "m",
                                'effective wavelength')]
    
    def _diminfo(self):

        nrows = len(self.data)
        return f"{super()._diminfo()}={nrows}W"

    def is_referred_to_by(self, other):
        return (not isinstance(other, _WavelengthHDU) and
                isinstance(other, _MustHaveWavelengthHDU) and
                self.get_insname() == other.get_insname())

    def __mod__(self, other):

        return False

    def rename(self, new_name):

        old_name = self.header['INSNAME']
        super().rename(new_name)
        
        container = self.get_container()
        if container is None:
            return

        for h in container.get_inspolHDUs():
            to_rename = h.data['INSNAME'] == oldname
            h.data['INSNAME'][to_rename][...] = new_name

    @classmethod
    def from_data(cls, *, insname, version=2, eff_wave, eff_band=0., 
            fits_keywords={}, **columns):

        """

        Build an OI_WAVELENGTH binary table from data

Arguments
---------

    insname (str)
        Name of this table
    version (int)
        OIFITS version if it cannot be deduced from context (optional)

    eff_wave (float, NWAVE)
        Effective wavelength
    eff_band (float, NWAVE)
        Bandwidth (optional, defaults to 0.)

    fits_keywords (dict)
        Additional FITS header keywords (optional)

Additional arguments
--------------------

Any additional keyword argument will be appended as a non-standard FITS 
column with its name prefixed with NS_ 

        """
        fits_keywords = dict(insname=insname, **fits_keywords)
        columns = dict(eff_wave=eff_wave, eff_band=eff_band, **columns)

        return super().from_data(fits_keywords=fits_keywords, **columns)

    def _trim_helper(self, *, wave_filter=lambda w: True, 
            target_filter=None, insname_filter=None, keep_ns_columns=False):

        keep = _np.vectorize(wave_filter)(self.EFF_WAVE)

        columns = {}
        standard_colnames = self._get_oi_colnames()

        for column in self.columns:

            colname = column.name
            if not keep_ns_columns and colname not in standard_colnames:
                continue
            columns[colname.lower()] = self.data[colname][keep]

        whdu = self._from_data(fits_keywords=self.header, **columns)
        
        return whdu
        

    def _bin_helper(self, R):

        min, max = _np.minimum, _np.maximum

        # convert to double precision (rounding errors in weights...)
        wave = _np.asarray(self.get_wave(shape='wavelength'), dtype=float)
        band = _np.asarray(self.get_band(shape='wavelength'), dtype=float)

        winf, wsup = wave - band / 2,  wave + band / 2
        wmin, wmax = winf[0], wsup[-1]
        wmed = (wmin + wmax) / 2
        R0 = wmed / ((wmax - wmin) / len(wave))

        if R0 > R and len(wave) > 1:

            nwave = int(_np.ceil(len(wave) * R / R0))

            dw = (wmax - wmin) / nwave
            w1 = winf[0] + dw / 2
            w2 = wsup[-1] - dw / 2
            new_wave = _np.linspace(w1, w2, nwave)
            new_winf = new_wave - dw / 2
            new_wsup = new_wave + dw / 2
            new_band = _np.full((nwave,), dw)

            weights = max(min(wsup, new_wsup[:,None])
                            - max(winf, new_winf[:,None]), 0) / band
            weights = weights.T
            
            oi_colnames = self._get_oi_colnames()
            fmt = f"{nwave}E"

            wave_col = _fits.Column('EFF_WAVE', fmt, 'm', array=new_wave) 
            band_col = _fits.Column('EFF_BAND', fmt, 'm', array=new_band) 
            other_cols = [c for c in self.columns if c.name not in oi_colnames]
            cols = [wave_col, band_col, *other_cols]

            new = self.from_columns(cols, header=self.header)

        else:
            
            new = self.copy()
            weights = None

        return new, weights



class WavelengthHDU1(
        _WavelengthHDU,
        _OITableHDU11,
     ):
    """

First revision of the OI_WAVELENGTH binary table, OIFITS v. 1

    """
    _COLUMNS = [('EFF_BAND', True, '1E', (), _u.is_strictpos, None, "m",
        'effective bandwidth')]


class WavelengthHDU2(
        _WavelengthHDU,
        _OITableHDU22,
      ):
    """

First revision of the OI_WAVELENGTH binary table, OIFITS v. 1

    """
    _COLUMNS = [('EFF_BAND', True, '1E', (), _u.is_pos, None, "m",
        'effective bandwidth')]
  
new_wavelength_hdu = _WavelengthHDU.from_data 
