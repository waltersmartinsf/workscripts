from astropy import units as u #units package
from astropy.coordinates import SkyCoord #coordinate system package
from astropy.coordinates import EarthLocation, AltAz #obtain a Earth Location Coordinate
from astropy.coordinates import Longitude, Latitude
from astropy.coordinates import Angle #work with angles
from astropy.coordinates import get_sun,get_moon #Sun and Moon Location on Sky
from astropy.time import Time
import ephem

def hourangle():
    '''
    Obtain the Hour Angle for a specific date, site and sky coordinates.
    ---
    INPUT:
    
    '''
    return HA 