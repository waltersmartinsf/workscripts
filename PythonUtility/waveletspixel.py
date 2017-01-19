'''
LIbrary to work with wavelets using:

--- Requirements ---
Python >2.7
Numpy
PyWavelets
Astropy
Matplotlib
---

Functions:

wavelets_fits: to run inside a directory, return the mean data and the details 
'''
def wavelets_fits(path,idx_im,wav,mode='sym',show=False):
    import pywt #python wavelets package
    import os
    from astropy.io import fits
    import numpy as np
    import glob
    import matplotlib.pyplot as plt
    
    images = glob.glob(path+idx_im+'*.fits')
    print('Working on those images = \n')
    media_images_wt = []
    details = []
    for i in images:
        print('Image = '+i+'\n')
        im,hdr = fits.getdata(i,header=True) #reading the fits image (data + header)
        im = np.array(im,dtype='Float64') #transform the data for a matrix
        tam = np.shape(im) #dimension of the matrix
        if show == True:
            plt.figure()
            plt.imshow(im,origin='lower', cmap=plt.cm.gray,vmin=np.mean(im)-np.std(im),vmax=np.mean(im)+np.std(im))
            plt.title('Original Image \n'+i,fontsize=10)
        (media_wt,(Hdet_wt,Vdet_wt,Ddet_wt)) = pywt.dwt2(im,wav,mode=mode)
        if show == True:
            plt.figure()
            plt.imshow(media_wt,origin='lower', cmap=plt.cm.gray,vmin=np.mean(media_wt)-np.std(media_wt),
                       vmax=np.mean(media_wt)+np.std(media_wt))
            plt.title('Wavelet Transform '+wav+'\n'+i,fontsize=10)
        media_images_wt.append(media_wt)
        details.append([(Hdet_wt,Vdet_wt,Ddet_wt)])
    return media_images_wt,details

def wavelets_fits_region(path,idx_im,wav,center,mode='sym',delta=7.,show=False):
    import pywt #python wavelets package
    import os
    from astropy.io import fits
    import numpy as np
    import glob
    import matplotlib.pyplot as plt
    
    images = glob.glob(path+idx_im+'*.fits')
    print('Working on those images = \n')
    images_data = []
    media_images_wt = []
    details = []
    for i in images:
        print('Image = '+i+'\n')
        im,hdr = fits.getdata(i,header=True) #reading the fits image (data + header)
        im = np.array(im,dtype='Float64') #transform the data for a matrix
        im = im[int(center[0])-delta:int(center[0])+delta,int(center[1])-delta:int(center[1])+delta]
        tam = np.shape(im) #dimension of the matrix
        images_data.append(im)
        if show == True:
            plt.figure()
            plt.imshow(im,origin='lower', cmap=plt.cm.gray,vmin=np.mean(im)-np.std(im),vmax=np.mean(im)+np.std(im))
            plt.title('Original Image \n'+i,fontsize=10)
        (media_wt,(Hdet_wt,Vdet_wt,Ddet_wt)) = pywt.dwt2(im,wav,mode=mode)
        if show == True:
            plt.figure()
            plt.imshow(media_wt,origin='lower', cmap=plt.cm.gray,vmin=np.mean(media_wt)-np.std(media_wt),
                       vmax=np.mean(media_wt)+np.std(media_wt))
            plt.title('Wavelet Transform '+wav+'\n'+i,fontsize=10)
        media_images_wt.append(media_wt)
        details.append([(Hdet_wt,Vdet_wt,Ddet_wt)])
    return media_images_wt,details,images_data