from numpy import isnan, cos, sin
import astropy.units as u
from astropy.coordinates import AltAz

def altaz_frame(location, *, obswl, refraction=False):
    """

Alt-az frame taking into account mean astmospheric refraction.

Arguments
---------

location (astropy.coordinates.EarthLocation)
    Earth location of site 

obswl (astropy.units.Quantity)
    Observing wavelength


Returns
-------

altaz (astropy.coordinates.AltAz)
    Alt-az frame

    """
    
    if not refraction:
        return AltAz(location=location)

    lat = location.lat.to_value("rad")
    h = location.height.to_value("m")

    # Mean temperature at sea level  
    # 
    # Rough approximation to Feulner et al. (2013). On the Origin of 
    # the Surface Air Temperature Difference between the Hemispheres 
    # in Earth's Present-Day Climate. Journal of Climate. 26. 7136-7150.
    # 10.1175/JCLI-D-12-00636.1. 

    T0 = 250 + 50 * cos(lat) # K

    # Mean temperature on site
    #
    # Neumann (1955). Latitudinal variation of tropospheric temperature
    # lapse rate. Arch. Met. Geoph. Biokl. A. 8, 351–353  
    # https://doi.org/10.1007/BF02247093 
    # But we use the 0.0065 value of the Standard US Atmosphere for
    # the low to moderate latitudes

    Lb = 6.5e-3 - 2e-3 * sin(lat) ** 4 # K/m
    T = (T0 - h * Lb - 273.15) # ⁰C 

    # Mean pressure at sea level
    #
    # Variation with latitude is negligible at this level of accuracy
    # (T, P vary easily by 5%).

    P0 = 101325 # Pa

    # Mean pressure on site
    #
    # Let's use the barometric formula for a real atmosphere 
    # (dT/dz < 0, hT = -dz/dT)  to estimate the pressure at 
    # the observatory. It will give an estimate for the 
    # atmospheric refraction needed to determine the apparent
    # position of the source
    
    h0 =  8438 * (T0 / 288.15) # m (isothermal scale height)
    hT = 44330 * (6.25e-3 / Lb) # m (temperature scale height)

    P = P0 * (1 - h/hT) ** (hT/h0) # Pa

    # Relative humidity
    H = 20 * u.percent # % (typical of observatories)

    frame = AltAz(location=location, 
                    pressure=P * u.Pa, 
                    temperature=T * u.deg_C, 
                    obswl=obswl,
                    relative_humidity=H * u.percent) 
    
    return frame

def apply_space_motion(coo, obstime, correct_motion=False):
    """

Update coordinates to current epoch

    """
    # set new epoch
    dt = obstime - coo.obstime
    coo.obstime = obstime
 
    # position correction due to motion (linear approximation)
    pmra = coo.pm_ra_cosdec
    pmde = coo.pm_dec
    rv = coo.radial_velocity
    d = coo.distance
    
    has_pm = ~(isnan(pmra) | isnan(pmde))
    has_rv = ~isnan(rv)

    coo.ra[has_pm] += pmra[has_pm] / cos(coo.dec[has_pm]) * dt[has_pm]
    coo.dec[has_pm] += pmde[has_pm] * dt[has_pm]
    coo.distance[has_rv] += rv[has_rv] * dt[has_rv]

    # motion correction
    # we probably don't need it for OIFITS...
    if correct_motion:
        raise RuntimeError('Not implemented')

    return coo
