'''
Created by Walter Martins Filho

This script have a set of functions to fitting a Gaussian or a pseudo-Voigt profile to a stellar PSF
and return the FWHM, for that night. Also, it have functions to normalize each PSF by its own standart
devition retrived by the fitting model.
'''
import numpy as np
import pandas as pd
import os
from glob import glob
import matplotlib.pylab as plt
from astropy.io import fits
import lmfit
import simulate_starfield as smpf
import matplotlib as mpl

def gaussian2d(x,y,x0,y0,amp,fwhm):
    return amp * np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

def star_psf(x,y,x0,y0,a,sigx,sigy,b):
    gaus = a * np.exp(-(x-x0)**2 / (2*sigx**2) ) * np.exp(-(y-y0)**2 / (2*sigy**2) ) + b
    return gaus

# def exposition(x0,y0,a,sigx,sigy,rot,b,CCDsize,show = False):
#     x0,y0 = int(x0),int(y0)
#     img = smpf.ccd(CCDsize)
#     star = smpf.psf(x0,y0,a,sigx,sigy,rot,b)
#     img.draw(star)
#     img.data = img.data + sky_bkg_model(CCDsize)

#     return img.data

def exposition(x0,y0,a,sigx,sigy,rot,b,CCDsize,show = False):
    def sky_bkg_model(CCDsize):
        field = 10+ np.zeros((CCDsize[0],CCDsize[1]))
        return np.random.poisson(field)
    img = smpf.ccd(CCDsize)
    star = smpf.psf(x0,y0,a,sigx,sigy,rot,b)
    img.draw(star)
    img.data = img.data + sky_bkg_model(CCDsize)

    pars_psf = smpf.fit_centroid(img.data,[x0,y0],box=25,psf_output=False)
    if show == True:
        print('centeroid bg=',pars_psf[-1])
    area = smpf.phot(pars_psf[0],pars_psf[1],img.data,r=15,debug=False,bgsub=True)
    if show == True:
        print(pars_psf)
        print('phot area=',area)
        print('psf area=',star.area())
    return img, star, pars_psf, area

def transit_model(t,Rp,aR,P,i,u1,u2,e,omega,tmid):

    import ctypes
    import requests

    n = len(t)

    # define C double array
    array_1d_double = np.ctypeslib.ndpointer(dtype=ctypes.c_double,ndim=1,flags=['C_CONTIGUOUS','aligned'])
    # load library
    lib_trans = np.ctypeslib.load_library('lib_transit.so','/Users/walter/github/workscripts/PythonUtility/lightcurveModel')
    # load function occultquadC from transit library
    occultquadC = lib_trans.occultquad
    occultquadC.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, array_1d_double]
    occultquadC.restype = None
    curve = np.zeros(n,dtype=ctypes.c_double)
    
    t = np.require(t,dtype=ctypes.c_double,requirements='C')
    occultquadC(t,Rp,aR,P,i,u1,u2,e,omega,tmid,n,curve) # old results
    return curve

def import_fits(image_path,Xc,Yc,delta):
    im0 = fits.getdata(image_path,header=False)
    im0 = np.array(im0,dtype='Float64')
    #cut the image in the close region the the XO-2b
    #img[columns:lines] 
    img = im0[int(Xc-delta):int(Xc+delta),int(Yc-delta):int(Yc+delta)]
    return img

def cut_image(img,Xc,Yc,delta):
    '''
    Cut a image in a numpy array object.
    '''
    img = img[int(Xc-delta):int(Xc+delta),int(Yc-delta):int(Yc+delta)]
    return img

def centerPSF(img):
    '''
    GOAL: 
    estimating the center of the PSF using the maximum values in the column and in the rows of a grid.

    This function only work with the grid have ONE, AND ONLY ONE, star in the grid.

    INPUT:
    
    img: 2D-array with double-float values

    OUTPUT:

    center: coordinates in pixel of that grid that have the maximum value.

    '''
    center = np.unravel_index(np.argmax(img),img.shape)
    return center

def gaussian(x, amp, cen, wid):
    '''
    gaussian PSF
    '''
    return amp * np.exp(-(x-cen)**2 /wid)

# PSF Fitting
def psf_model(rows,counts,model='gaussian', amp = 5000, cen = 80, wid = 1., W = 1., B=1., show=False):
    '''
    Fitting a PSF model to a column choosing between the Gaussian or pseudo-Voigt profile.
    '''
    def gaussian(x, amp, cen, wid):
        '''
        gaussian PSF
        '''
        return amp * np.exp(-(x-cen)**2 /wid)
    
    def pvoigt(x, amp, cen, sigma, W, B):
        '''
        pseudo-Voigt profile
        PSF model based on Kyle's paper
        '''
        return (1-W) * amp * np.exp(-(x - cen)**2/(2.*sigma)) + 1. * W  * (amp*sigma**2)/((x-cen)**2 +sigma**2) + B
    if model == 'gaussian':
        gmodel = lmfit.Model(gaussian)
        result = gmodel.fit(counts, x=rows, amp=amp, cen=cen, wid=wid)
        
    if model == 'pVoigt':
        gmodel = lmfit.Model(pvoigt)
        result = gmodel.fit(counts, x=rows, amp=amp, cen=cen, sigma = wid/2., W = 1., B = 1.) 

    if show == True:
        print(result.fit_report())
    return result


def psf_model_plot(rows, counts, result, save_dir, save_name, show = False):
    fwhm = float(2.355 * np.sqrt(result.params['wid'].value/2.))
    plt.figure()
    # plt.plot(rows, counts,'bo',color='red',label='psf')
    # plt.plot(rows, counts,'bo',color='red',label='psf')
#     plt.plot(rows, result.init_fit, 'k--')
    plt.plot(rows, result.best_fit, 'b-',label='model')
    # plt.vlines(result.params['cen'].value+result.params['wid'].value/2.,counts.min()-100,counts.max()+100)
#     plt.vlines(int(result.params['cen'].value),counts.min(),counts[np.ceil(result.params['cen'].value)],
#                linestyles='dashed')
    plt.hlines(counts.min(),rows.min(),rows.max(),linestyles='dashed')
    plt.hlines(counts.max()/2.,result.params['cen'].value-fwhm/2.,result.params['cen'].value+fwhm/2.,color='red')
    # plt.text(result.params['cen'].value,counts.max()/2.,'FWHM = '+str(round(fwhm,4)))
    # plt.title('FWHM = '+str(round(fwhm,4)))
    plt.title('FWHM = '+str(9.420))
    plt.legend()
    plt.xlabel('column')
    plt.ylabel('counts')
    plt.savefig(save_dir+save_name)
    if show == True:
        plt.show()
    plt.close()

def flux_aperture(image,Xc,Yc,sigx):
    from photutils import CircularAperture
    from photutils import aperture_photometry
    from photutils import CircularAnnulus

    #local where we will apply th photometry
    apertures = CircularAperture((Xc,Yc), r=sigx)
    annulus_apertures = CircularAnnulus((Xc,Yc), r_in=sigx+2, r_out=sigx+4)
    #perform photometry
    apers = [apertures, annulus_apertures]
    phot_table = aperture_photometry(image, apers)
    # phot_table = aperture_photometry(image- bkg, apertures)
    #subtracted the sky background
    bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
    bkg_sum = bkg_mean * apertures.area()
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    
    return final_sum

def normalize_psf(img,fwhm,center,amp):
    sigma = fwhm/2.355
    rawflux = flux_aperture(img,int(center[0]),int(center[1]),sigma)
    # img_norm = (img-np.mean(img))/np.std(img)
    img_norm = abs(img/rawflux)
    fit_model = psf_model(np.arange(0,img.shape[0],1),img_norm[int(center[0]),:],model='gaussian', amp = 0.1, cen = center[1], wid = 2*sigma**2, show=False)
    # img_norm = (img/rawflux)*gaussian(x=np.arange(0,img.shape[0],1),amp=amp/rawflux,cen=center[0],wid=2.*sigma**2)
    fwhm = 2.355 * np.sqrt(fit_model.params['wid']/2.)
    return img_norm, fwhm

def convolute_psf(img0,img):
    '''
    This function convolute the first psf with the second one, make both have the same spatial scale.
    '''
    # Here, we will turn the counts in the image to small numbers:

    # img = img - img[centerXY[1],centerXY[0]]/2.
    # img0 = img0 - img0[centerXY[1],centerXY[0]]/2.
    # convolute both images
    img = img * img0
    #img = absimg - img.mean()
    return img

def ica_to_data(dataframe,path_to_dir):
    '''
    Apply Independent Component Analysis to a specific dataframe. 
    '''
    import scipy.io as sio #read matlab binaries
    import ica_utility
    from sklearn.decomposition import PCA
    import scipy
    import subprocess
    
    if os.path.isfile(path_to_dir+'data_to_ica.csv') == False:
        dataframe.to_csv(path_to_dir+'data_to_ica.csv')

    if os.path.isfile(path_to_dir+'ica_macosx.m') == False:
        os.system('cp /Users/walter/github/workscripts/ICAsource/ica_macosx.m '+path_to_dir)
    subprocess.call(["/Applications/MATLAB_R2016b.app/bin/matlab",'-nodisplay','-nosplash','-nodesktop','-r',"try, run('ica_macosx.m'), catch, exit, end, exit"])

    matlab_binaries = glob('*.mat')
    ica_data = sio.loadmat(matlab_binaries[0])
    print(ica_data.keys())
    signal = pd.DataFrame(ica_data['result']).T
    
    #Saving ISR information ********************************************************************
    ISR_xticks = np.arange(0,signal.shape[0],1,dtype='int')
    ISR_yticks = np.arange(0,signal.shape[1],1,dtype='int')
    Norm = (ica_data['ISRwasobi'].mean() + ica_data['ISRefica'].mean()) / 2.

    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True,figsize=(20,7))
    
    # plt.grid(color='white')
    algorithm = ['ISRwasobi','ISRefica']
    algorithm_name = ['ISR WASOBI','ISR EFICA']
    for i in range(len(axes.flat)):
        print(i)
        print(algorithm[i])
        im = axes.flat[i].imshow(ica_data[algorithm[i]]/Norm, vmin= 0, vmax= 5,aspect='auto')
        axes.flat[i].set_xticks(ISR_yticks)
        axes.flat[i].set_yticks(ISR_yticks)
        axes.flat[i].set_title(algorithm_name[i])

    cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
    plt.colorbar(im, cax=cax, **kw).set_label('Normalized at '+str(round(Norm,4)))
    plt.savefig(path_to_dir+'ISR_comparing_.png')
    
    #Check the lightcurve transit **************************************************************
    pca = PCA(n_components=len(signal.columns))
    component_id = ica_utility.pearson_ica_test(signal,dataframe,'./')

    if (signal[component_id][int(signal.index.values.mean())] > signal[component_id].mean()) == True:
        signal[component_id] = ica_utility.rotate_axis(signal[component_id])
    
    transit_signal = signal[component_id]

    return signal, ica_data, transit_signal

def pixelICA(dataframe,save_path,matlab_path="/Applications/MATLAB_R2016b.app/bin/matlab"):
    import scipy.io as sio #read matlab binaries
    import ica_utility
    from sklearn.decomposition import PCA
    import scipy
    import subprocess
    # import sourcesfits
    # from photutils import aperture_photometry, Background, CircularAperture
    # from astropy.stats import sigma_clipped_stats #find psf sources
    # from photutils import daofind
    workdir = os.getcwd()
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    os.chdir(save_path)
    dataframe.to_csv('data_to_ica.csv')
    if os.path.isfile('ica_macosx.m') == False:
        os.system('cp /Users/walter/github/workscripts/ICAsource/ica_macosx.m ./')
    subprocess.call(["/Applications/MATLAB_R2016b.app/bin/matlab",'-nodisplay','-nosplash','-nodesktop','-r',"try, run('ica_macosx.m'), catch, exit, end, exit"])
    matlab_binaries = glob('*.mat')
    ica_results = sio.loadmat(matlab_binaries[0])
    
    plt.figure()
    plt.title('ICA Results')
    plt.imshow(ica_results['result'].T,origin='lower')
    plt.colorbar()
    plt.savefig(save_path+'ICA_results.png')

    #Check the lightcurve transit **************************************************************
    ica_result_df = pd.DataFrame(ica_results['result'])

    #creating a transit model to compare
    time = pd.read_csv('/Users/walter/OneDrive/work/ICA_Analysis/Simulations/time_results.csv')
    time = (time.HJD - time.HJD.mean())
    model = transit_model(time.values,0.1037,8.2332,2.615838,88.7,0.81461289, 0.02619382,0.045,0,0)
    model_transit = pd.DataFrame()
    contador = 0
    while contador < len(dataframe.columns):
        xx = model + np.random.normal(0,0.0001,size=len(model))
        model_transit = pd.concat([model_transit,pd.DataFrame(xx)],axis=1)
        contador += 1

    # component_id = ica_utility.pearson_ica_test(ica_result_df.T,model_transit.T,'./')
    ica_result_df = ica_result_df.T #turn each column in a a flux
    dataframe = dataframe.T
    component_id = ica_utility.pearson_ica_test(ica_result_df,dataframe,'./')
    if (ica_result_df[component_id][int(ica_result_df.index.values.mean())] > ica_result_df[component_id].mean()) == True:
        ica_result_df[component_id] = ica_utility.rotate_axis(ica_result_df[component_id])
    transit_signal = ica_result_df[component_id]
    

    plt.figure()
    plt.plot(transit_signal)
    plt.savefig(save_path+'transit_from_pixelICA.png')
    os.chdir(workdir)
    return ica_results,transit_signal

def ica_df(img,delta):
    '''
    This function cut the image and transform a matrix object in a row object to be applied in an ICA routine
    '''
    centerXY = centerPSF(img)
    img = cut_image(img,centerXY[0],centerXY[1],delta)
    rows,cols = img.shape[0],img.shape[1]
    img = img.flatten()
    img = pd.DataFrame(img)
    return img, rows, cols

def fwhm_xo2b(length,data_night):
    '''
    Simulate the FWHM from the XO-2b's data set using a spline function combined with a normal noise.
    '''
    from scipy.interpolate import interp1d
    # import scipy as sp

    fwhm = pd.read_csv(data_night,header=None)
    y = fwhm[0].values
    x = np.linspace(0,1,len(y))
    new_length = length
    new_x = np.linspace(x.min(), x.max(), new_length)
    #new FWHM for the night:
    new_y = interp1d(x, y, kind='cubic')(new_x)+np.random.normal(0.1,0.01,size=new_length)
    #new sigma:
    new_sigma = (new_y/2.355)
    return new_sigma

def sky_bkg_model(CCDsize):
    field = 10+ np.zeros((CCDsize[0],CCDsize[1]))
    return np.random.poisson(field)

# def seeing_model(size):
#     t = np.linspace(0,1,size)
#     return np.exp(t)

def mixing_signal(unmixed_signal,size,noise='sinoidal',scale=10):
    import scipy.stats as stats

    step = np.linspace(0,1,size)

    if noise == 'linear':
        signal2 = (-0.0001*step+0.0001)*scale
    # S = np.c_[unmixed_signal, signal2]
    # S -= S.mean(axis=0)
    # S /= S.std(axis=0)
    # A = np.array([[0.5, 0.5], [0.5, 0.5]])  # Mixing matrix Test
    # X = np.dot(S, A.T)
    # X = pd.DataFrame(X)
    # return X[0].values
    if noise == 'sinoidal':
        signal2 = scale*0.001*np.sin(2*np.pi*step/0.5)
    
    if noise == 'normal':
        signal2 = np.random.normal(0,0.0001*scale,size)
    
    if noise == 'exponential':
        signal2 = scale*0.001*np.exp(-step)
    
    if noise == 'gaussin_drift':
        signal2 = scale*0.0001*stats.norm(0.3,0.1).pdf(step)

    if noise == 'twomodes':
        signal2 = scale*0.0001*(stats.norm(0.3,0.1).pdf(step) + stats.norm(0.7,0.1).pdf(step))
        
    return unmixed_signal + signal2

def psfdrift(x0,y0,size):
    '''
    Create the drift movement for a point spread function in the CCD during the night, for N-size expositions.
    ---
    x0, y0: initial positions in the CCD
    size: number of expositions along the night
    '''
    # x_drift,y_drift = np.random.normal(size=size),np.random.normal(size=size)
    x_drift,y_drift = np.random.choice([0,1,-1],size=size),np.random.choice([0,1,-1],size=size)
    xx,yy = np.zeros(size), np.zeros(size)
    for i in range(size):
        if i == 0:
            xx[i],yy[i] = x0 + x_drift[0], y0 + y_drift[0]
        else:
            xx[i],yy[i] = xx[i-1]+x_drift[i], yy[i-1]+y_drift[i]
    return xx, yy

def fitspline(x,size):
    from scipy.interpolate import interp1d
    y = x
    x = np.linspace(0,1,len(y))
    new_length = size
    new_x = np.linspace(x.min(), x.max(), new_length)
    new_y = interp1d(x, y, kind='cubic')(new_x)
    return new_y

def fwhm_noise(sigma_x, size, scale=1, kind='normal'):
    '''
    Create a noise fwhm in a constante input value over time.
    '''
    if kind == 'normal':
        noise = scale*np.random.normal(0,sigma_x,size)
    return noise