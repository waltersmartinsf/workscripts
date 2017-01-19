import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame
# import os
# import glob

from astropy.io import fits
# from astropy.wcs import WCS
# from astropy.coordinates import Angle
# from astropy import units as u
# from astropy.coordinates import SkyCoord

from photutils import datasets
from astropy.stats import sigma_clipped_stats
from photutils import daofind

# from photutils import CircularAperture
# from astropy.visualization import SqrtStretch
# from astropy.visualization.mpl_normalize import ImageNormalize

def find_sources_fits(image,fwhm,show=False):
	if show == True:
		print('Obtain the point source locations at the FITS image = \n')
		print(image)
		print('\n')
	im,hdr = fits.getdata(image,header=True) #reading the fits image (data + header)
	im = np.array(im,dtype='Float64') #transform the data for a matrix
	tam = np.shape(im) #dimension of the matrix
	mean, median, std = sigma_clipped_stats(im, sigma=fwhm, iters=5)
	if show == True:
		print('Mean, Median, STD = \n')
		print(mean, median, std)
		print('\n')
	sources = daofind(im - median,fwhm=fwhm, threshold=5.*std)
	if show == True:
		print('Sources found! \n')
		print(sources)
		print('\n')
		plt.figure()
		plt.imshow(im,origin='lower', cmap=plt.cm.gray,vmin=np.mean(im)-np.std(im),vmax=np.mean(im)+np.std(im))
		plt.colorbar()
		plt.scatter(sources['xcentroid'],sources['ycentroid'],color='red')
	#Transform astropy table in a pandas dataframe
	data = []
	for i in sources.keys():
		data.append(sources[i])
	data = DataFrame(data).T
	data.columns = sources.keys()
	return data 

def find_sources_array(image,fwhm,show=False):
	if show == True:
		print('Obtain the point source locations at the FITS image = \n')
		print(image)
		print('\n')
	# im,hdr = fits.getdata(image,header=True) #reading the fits image (data + header)
	# im = np.array(im,dtype='Float64') #transform the data for a matrix
	tam = np.shape(image) #dimension of the matrix
	mean, median, std = sigma_clipped_stats(image, sigma=fwhm, iters=5)
	if show == True:
		print('Mean, Median, STD = \n')
		print(mean, median, std)
		print('\n')
	sources = daofind(image - median,fwhm=fwhm, threshold=5.*std)
	if show == True:
		print('Sources found! \n')
		print(sources)
		print('\n')
		plt.figure()
		plt.imshow(image,origin='lower', cmap=plt.cm.gray,vmin=np.mean(image)-np.std(image),
			vmax=np.mean(image)+np.std(image))
		plt.colorbar()
		plt.scatter(sources['xcentroid'],sources['ycentroid'],color='red')
	#Transform astropy table in a pandas dataframe
	data = []
	for i in sources.keys():
		data.append(sources[i])
	data = DataFrame(data).T
	data.columns = sources.keys()
	return data 