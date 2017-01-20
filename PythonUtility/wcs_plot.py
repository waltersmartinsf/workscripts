"""

Goal: wcs projection

Project the World Coordinate System on a plot image
---
Author: Walter S. Martins-Filho
email:  walter@on.br
        waltersmartinsf@gmail.com
---

"""

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

def wcs_projection(fits_path):
    im, hdr = fits.getdata(fits_path, header=True)
    im = im = np.array(im,dtype='Float64') #transform the data for a matrix
    tam = np.shape(im) #dimension of the matrix
    wcs = WCS(hdr)
    fig = plt.figure()
    fig.add_subplot(111, projection=wcs)
    plt.imshow(im, origin='lower', cmap=plt.cm.viridis)
    plt.xlabel('RA')
    plt.ylabel('Dec')
