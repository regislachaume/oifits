from .table import _OITableHDU
    
class _Referenced(_OITableHDU):
    """OI tables that are explicitely referenced to by other ones (OI_ARRAY, 
OI_WAVELENGTH, OI_CORR)."""
    
    def get_referrers(self):
        """Get all tables referring to this table."""
        
        container = self.get_container()
        if not container:
            return []

        key = self._REFERENCE_KEY
        val = self.header[key]

        refs = [h for h in container 
                    if h.header.get(key, '') == val and h is not self ] 

        return refs

    def rename(self, new_name):
        """Rename the table and updates all tables referring to it"""
        
        refkey = self._REFERENCE_KEY
        container = self.get_container()
        
        if container:

            for hdu in container:
                if hdu.header.get(refkey, None) == new_name:
                    raise RuntimeError("{refkey} = {new_name} already in use")

            for hdu in self.get_referrers():
                hdu.header[refkey] = new_name

        self.header[refkey] = new_name
