#Created by Kyle Pearson

from astropy.io import fits
import numpy as np

def merge_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def ra2deg(rasta):
    '''
    conver string '17:57:57.69' into degrees
    '''
    slist = rasta.split(":")
    hour = float(slist[0])
    minute = float(slist[1])
    second = float(slist[2])
    deg = 15*hour + minute*(15./60.) + second*15./(60*60)
    return deg

def dec2deg(decstr):
    slist = decstr.split(":")
    dec = float(slist[0])
    minute = float(slist[1])
    second = float(slist[2])
    deg = dec + minute*(1/60.) + second*1./(60*60)
    return deg


class image(object):
    def __init__(self,size):

        # load data from array
        if isinstance(size,np.ndarray):
            self.data = np.copy(size)
        else:
            self.data = np.zeros(size)

        self.header = {
            "OBSERVER" : "Observer",
            "INST"     : "NESSI",
            "TELESCOP" : "P200in",
            "FILENAME" : "default",
            "IMGTYPE"  : "fits",
            "RA"       : "22:03:10.77207", #update
            "DEC"      : "+18:53:03.5430", #update

            #gnomonic projection http://docs.astropy.org/en/stable/wcs/#supported-projections
            # WCS keywords below
            "CTYPE1"   : "RA---TAN",
            "CTYPE2"   : "DEC--TAN",    #gnomonic - same as HUBBLE WFC
            "CRPIX1"   : 1024,         # ref point pixel x
            "CRPIX2"   : 1024,         # ref point pixel y
            "CDELT1"   : 0.5*0.000119444444,# deg per pixel x (0.02)
            "CDELT2"   : 0.5*0.000119444444,# deg per pixel y
            "CRVAL1"   : 0, # right ascension of ref pixel (deg)
            "CRVAL2"   : 0,  # get from telescope
            "CROTA2"   : 0 # Image rotation value # add PA from telescope?
            }
        self.header['CRVAL1'] = ra2deg(self.header['RA'])
        self.header['CRVAL2'] = dec2deg(self.header['DEC'])

    def save(self,name='test.fits'):
        header = fits.header.Header(self.header)
        hdu = fits.PrimaryHDU(self.data,header=header)
        hdu.writeto(name,overwrite=True)

if __name__ == "__main__":

    # open 61" image
    # pass data to image class
    # merge headers using function above
    # update header key words with valid WCS system (http://docs.astropy.org/en/stable/wcs/#building-a-wcs-structure-programmatically)

    # THINGS TO UPDATE
    # CRPIX1 = 646.8 # centroid of WASP-69 on 61" image, this is valid pixel
    # CRPIX2 = 591.7
    # CRVAL1  = ra2deg('21:00:06.19')
    # CRVAL2 = dec2deg('-5:05:40.1')
    # CROTA2 - CHECK image header for position angle

    # CDELT1 = PIXSCAL1*3 # 3 is for 3x3 binning from image header
    # CDELT2 = PIXSCAL2*3 # 3 is for 3x3 binning from image header

    # after a new image is saved, open it in ds9
    # verify position (RA,DEC) of other stars with simbad, use ds9 to get RA and DEC of each pixel
    # check that it can be opened with the SOAR software
