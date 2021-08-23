from .table import _OITableHDU

__all__ = []
    
class _Referenced(_OITableHDU):
    """OI tables that are explicitely referenced to by other ones (OI_ARRAY, 
OI_WAVELENGTH, OI_CORR)."""
    
    def _get_referrers(self, cls):
        """Get all data tables referring to this table"""
    
        if (container := self.get_container()) is None:
            return []
        
        key = self._REFERENCE_KEY
        val = self.header[key]

        def has_same_key(hdu):
            return hdu is not self and hdu.header.get(key, None) == val

        hdus = container.get_HDUs(cls, filter=has_same_key)

        return hdus

    def get_referrers(self):
        """Get all tables referring to this table."""
        
        from .inspol import _InspolHDU      
 
        return self._get_referrers((_OITableHDU, _InspolHDU)) 

    def get_data_referrers(self):
        """Get all data tables (OI_FLUX, OI_VIS, OI_VIS2, OI_T3) referring
        to this table"""
        
        # necessarily here to avoid circular def.
        from .data import _DataHDU

        return self._get_referrers(_DataHDU)

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
