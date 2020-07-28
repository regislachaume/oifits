from astropy.io import fits

class PrimaryHDU(fits.PrimaryHDU):
    def __repr__(self):

        if self.data:
            return f"<PrimaryHDU ({'×'.join(self.data.shape)} image)>"
        return '<PrimaryHDU>'
