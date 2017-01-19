"""
script to apply the independent component analysys on fits images using Astropy and
Numpy packages to work with matrices, and Pandas package to work with dataframe

Author: Walter S. Martins-Filho
email:  walter@on.br
        waltersmartinsf@gmail.com
___

This code uses a MATLAB routines created by Ingo Waldmann

___

References:

Waldmann (2012) - OF “COCKTAIL PARTIES” AND EXOPLANETS

"""
#*******************************************************************************
#*******************************************************************************

# Packages useful

import os #control the OS commands like create a folder or move a file
import sys #package to help the OS-package
import subprocess #package that call another program, in this case the MATLAB
import yaml #input data without any trouble
from glob import glob #create a list of files with the same prefix and/or sufix
import time #control and save information about time

import numpy as np #package to work with vectors and matrices
import pandas as pd #package to work with dataframes and big data

import matplotlib as mpl #package to plot
import matplotlib.pyplot as plt #plot library "like" MATLAB
from mpl_toolkits.mplot3d import Axes3D #3D-projection

from astropy.io import fits #import FITS files to work in a Python ambient.

#*******************************************************************************
#*******************************************************************************

def find_psf(image,fwhm,show=False):
	'''
	Obtain the position of the centroids in a image matrix
	---
	Input:

	- image:	str		image path including the name
	-
	'''
	#importing useful packages from astropy
	from photutils import datasets
	from astropy.stats import sigma_clipped_stats
	from photutils import DAOStarFinder #daofind

	im,hdr = fits.getdata(image,header=True) #reading the fits image (data + header)
	im = np.array(im,dtype='Float64') #transform the data for a matrix
	tam = np.shape(im) #dimension of the matrix
	mean, median, std = sigma_clipped_stats(im, sigma=fwhm, iters=5)

	# sources = daofind(im - median,fwhm=fwhm, threshold=5.*std)
	result = DAOStarFinder(threshold=median,fwhm=3.5)
	sources = result.find_stars(im)

	if show == True:
		plt.figure()
		plt.imshow(im,origin='lower', cmap=plt.cm.gray,vmin=np.mean(im)-np.std(im),vmax=np.mean(im)+np.std(im))
		plt.colorbar()
		plt.scatter(sources['xcentroid'],sources['ycentroid'],color='red')
		plt.savefig('psf_sources_.png')

	# sources = sources.to_
	return sources
#*******************************************************************************
#*******************************************************************************

def import_fits(index_image,data_path,workdir,show=True,save=True,sources=False,cut_psf=False,cut_center=(0,0),cut_delta=4):
	'''
	Import FITS images to a matrix 64bits variable
	---
	Input:

	- show:				boolean value to show to stop the code and show previous results.
	- save:				boolean value to save previous results

    ---
    Output:

    Return the list of images and a list of matrices where each element is an fits image in matrix format.
	'''

	matlab_path,data_path,workdir,index_image,location_package = info_yaml(yaml_path)
	if show == True:
		print(matlab_path,data_path,workdir,index_image)

	os.chdir(data_path)
	images_list = glob(index_image+'*.fits')

	#check the FoV of the first image and go to the work diretory
	if show == True:
		im, hdr = fits.getdata(images_list[0],header=True)
		im = np.array(im,dtype='Float64') #transform the data for a matrix
		tam = np.shape(im) #dimension of the matrix
		plt.figure()
		plt.imshow(im,origin='lower', cmap=plt.cm.gray,vmin=np.mean(im)-np.std(im),vmax=np.mean(im)+np.std(im))
		plt.colorbar()
		if save == True:
			os.chdir(workdir)
			plt.savefig('FoV_'+index_image+'_.png')
		else:
			os.chdir(workdir)
		plt.show()
	else:
		os.chdir(workdir)

	#cut PSF from the original images_list
	if sources == True:
		os.chdir(data_path)
		psf_sources = find_psf(images_list[0],fwhm=3.5,show=show)
		os.chdir(workdir)

	#cut the image in the region of the psf
	if cut_psf == True:
		os.chdir(data_path)
		Xc,Yc = cut_center
		images_cuts = []
		for i in images_list:
			im, hdr = fits.getdata(i,header=True)
			im = np.array(im,dtype='Float64') #transform the data for a matrix
			# cut_delta = int(cut_delta)
			im0 = im[int(Yc)-cut_delta:int(Yc)+cut_delta,int(Xc)-cut_delta:int(Xc)+cut_delta]
			images_cuts.append(im0)
		os.chdir(workdir)
		if save == True:
			for i in range(len(images_cuts)):
				plt.figure()
				plt.imshow(images_cuts[i],origin='lower', cmap=plt.cm.gray,vmin=np.mean(images_cuts[i])-np.std(images_cuts[i]),vmax=np.mean(images_cuts[i])+np.std(images_cuts[i]))
				plt.colorbar()
				plt.savefig(images_list[i]+'_cut_.png')
				plt.close()
		return images_list, images_cuts
    #if it is not to cut the image
	else:
        os.chdir(data_path)
		Xc,Yc = cut_center
		images_data = []
		for i in images_list:
			im, hdr = fits.getdata(i,header=True)
			im = np.array(im,dtype='Float64') #transform the data for a matrix
			# cut_delta = int(cut_delta)
			im0 = im[int(Yc)-cut_delta:int(Yc)+cut_delta,int(Xc)-cut_delta:int(Xc)+cut_delta]
			images_data.append(im0)
		return images_list, images_data


#*******************************************************************************
#*******************************************************************************
