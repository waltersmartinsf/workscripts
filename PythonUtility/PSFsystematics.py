import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import os
import glob
import time as Time
from astropy.io import fits
import scipy.stats as stats
from photutils import CircularAperture, aperture_photometry
import simulate_starfield as smpf #simulate a stellar PSF
import fwhmfit #fitting and change the FWHM of a point source
import update #update code bar

def ica_to_data(dataframe,save_dir,matlab_ica_script):
    import scipy.io as sio #read matlab binaries
    import scipy
    import subprocess
    
    workdir = os.getcwd()
    
    if os.path.isfile(save_dir+'data_to_ica.csv') == False:
        dataframe.to_csv(save_dir+'data_to_ica.csv')
        
    if os.path.isfile(save_dir+'ica_macosx.m') == False:
        os.system('cp '+matlab_ica_script+' '+save_dir)
        
    os.chdir(save_dir)
    subprocess.call(["/Applications/MATLAB_R2016b.app/bin/matlab",'-nodisplay','-nosplash',
                     '-nodesktop','-r',"try, run('ica_macosx.m'), catch, exit, end, exit"])
    
    matlab_binaries = glob.glob(save_dir+'*.mat')
    ica_data = sio.loadmat(matlab_binaries[0])
    os.chdir(workdir)
    
    return ica_data

def ica_analysis(ica_data,save_data,dataframe):
    import ica_utility
    from sklearn.decomposition import PCA
    
    ISR_xticks = np.arange(0,len(ica_data['result']),1)
    
    ica_data['ISRwasobi'] = np.real(ica_data['ISRwasobi'])
    ica_data['ISRefica'] = np.real(ica_data['ISRefica'])
    
    Norm = (ica_data['ISRwasobi'].mean() + ica_data['ISRefica'].mean()) / 2.

    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True,figsize=(20,7))
    # plt.grid(color='white')
    algorithm = ['ISRwasobi','ISRefica']
    algorithm_name = ['ISR WASOBI','ISR EFICA']
    for i in range(len(axes.flat)):
        print(i)
        print(algorithm[i])
        im = axes.flat[i].imshow(ica_data[algorithm[i]]/Norm, vmin= 0, vmax= 5,aspect='auto',cmap=plt.cm.coolwarm)
        axes.flat[i].set_xticks(ISR_xticks)
        axes.flat[i].set_yticks(ISR_xticks)
        axes.flat[i].set_title(algorithm_name[i])
        #including labels on the plot
        for j in range(ica_data[algorithm[i]].shape[0]):
            for k in range(ica_data[algorithm[i]].shape[1]):
                #k: linha
                #j: coluna
                #print(k,j,ica_data[algorithm[i]][j][k]/Norm)
                axes.flat[i].text(k,j,str(round(ica_data[algorithm[i]][j][k]/Norm,4)),
                              horizontalalignment='center',verticalalignment='center',color='white')

    cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
    plt.colorbar(im, cax=cax, **kw).set_label('Normalized at '+str(round(Norm,4)))
    plt.savefig(save_data+'ISR_comparing_.png')
    plt.savefig(save_data+'ISR_comparing_.pdf')
    plt.close()
    
    #Check the lightcurve transit **************************************************************
    signal = pd.DataFrame(ica_data['result']).T
    
    pca = PCA(n_components=len(signal.columns))
    component_id = ica_utility.pearson_ica_test(signal,dataframe,'./')

    if (signal[component_id][int(signal.index.values.mean())] > signal[component_id].mean()) == True:
        signal[component_id] = ica_utility.rotate_axis(signal[component_id])
    
    transit_signal = signal[component_id]
    
    return transit_signal
    
def systematic_amplitude(save_dir,time_result_file,planet_info,obs_points,scale=10.,kind='normal',
               save_plot=True,close_plot=True,simple_noise=False,):
    '''
    This function create the amplitude of PSF (normalized at the maximum) for each observation 
    for a simulated night. The amplitude will be convoluated with the maximum number of photons, 
    which will return the total flux over time. This assumes that the FWHM is constant over time.
    '''
    
    Rp,aR,P,i,u1,u2,e,omega,tmid = planet_info
    
    time = pd.read_csv(time_result_file)
    time = (time.HJD - time.HJD.mean())
    time = fwhmfit.fitspline(time.values,size)
#     print(time)
    #fwhmfit.transit_model(time,Rp,aR,P,i,u1,u2,e,omega,tmid)
    depth_model = fwhmfit.transit_model(time,Rp,aR,P,i,u1,u2,e,omega,tmid)
#     depth_model = fwhmfit.mixing_signal(depth_model,obs_points,noise=kind,scale=scale)

    #creating systematic noise -- ['linear','exp','sin','normal','twopeaks']
    if kind == 'normal':
        sys_noise = scale* Rp * np.random.normal(0,0.1,obs_points)
    if kind == 'sin':
        T = (time[-1]-time[0])/2.
        sys_noise = scale * Rp * np.sin(2*np.pi*time/T)
    if kind == 'linear':
        sys_noise = scale * Rp * time
    if kind == 'onepeak':
        step = np.linspace(0,1,obs_points)
        sys_noise = scale * Rp * stats.norm(0.3,0.1).pdf(step)
    if kind == 'twopeaks':
        step = np.linspace(0,1,obs_points)
        sys_noise = scale * Rp * (stats.norm(0.3,0.1).pdf(step) + stats.norm(0.7,0.1).pdf(step))
    if kind == 'exp':
        step = np.linspace(0,1,obs_points)
        sys_noise = scale* Rp* np.exp(-step)
        
    if save_plot == True:
        plt.figure(figsize=(9,9))
        plt.plot(time,sys_noise)
        plt.title('Systematic Noise: '+kind+' with scale '+str(scale)+r' $\times$ depth')
        plt.xlabel('Phase')
        plt.savefig(save_dir+kind+'_scale_'+str(scale)+'_'+'sys_noise.png')
        if close_plot == True:
            plt.close()
        else:
            plt.show()
            
    #add some normal "simple" noise
    if simple_noise == True:
        depth_model = depth_model + sys_noise + np.random.normal(0,0.001,obs_points)
    else:
        depth_model = depth_model + sys_noise
        
    if save_plot == True:
        plt.figure(figsize=(9,9))
        plt.plot(time, depth_model)
        plt.xlabel('Phase')
        plt.ylabel('Depth [Rp/Rs]')
        plt.title('systematic noise scale = '+str(scale)+r' $\times$ depth')
        plt.savefig(save_dir+kind+'_scale_'+str(scale)+'_'+'depth_mixing.png')
        if close_plot == True:
            plt.close()
        else:
            plt.show()
    
    #creating the amplitude
    amplitude = a*depth_model

    plt.figure(figsize=(9,9))
    plt.plot(time,amplitude)
    plt.xlabel('Phase')
    plt.ylabel('# counts')
    plt.savefig(save_dir+kind+'_scale_'+str(scale)+'_'+'amplitude.png')
    if close_plot == True:
        plt.close()
    else:
        plt.show()
    
    return depth_model,amplitude

def drift_model(save_dir,x0,y0,size):
    xx,yy = x0 + np.zeros(size), y0 + np.zeros(size)

    plt.figure(figsize=(9,9))
    plt.plot(xx,yy)
    plt.scatter(x0,y0,color='red')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(save_dir+'drift_.png')
    plt.close()
    return xx,yy

# drift_model(save_dir_orig,x0,y0,len(amplitude))

def circular_apperture(save_dir,img):
     #creating and alocate space for some variables
    rawflux,signal_ = np.zeros(size), np.zeros(size)
            
    df_rawflux,df_signal = pd.DataFrame(), pd.DataFrame()
    ica_norm_df, img_lines, img_columns = pd.DataFrame(), np.zeros(size), np.zeros(size)
    midpointX, midpointY = np.zeros(size),np.zeros(size)
    amplitude_fit = np.zeros(size)
            
    #width of the Gaussian PSF model
    wid = np.zeros(size)
    #fwhm retrived from the model
    fwhm = np.zeros(size)
    fwhm_norm = np.zeros(size)
    
    centerXY = fwhmfit.centerPSF(img)
    
    fit_result =  fwhmfit.psf_model(np.arange(0,img.shape[0],1),
                                        img[int(centerXY[0]),:],model='gaussian', 
                                        amp=img[int(centerXY[0]),:][img[int(centerXY[0]),:].argmax()]/2.,
                                        cen = centerXY[1], wid = 1.,show=False)
    fwhmfit.psf_model_plot(np.arange(0,img.shape[0],1),img[int(centerXY[0]),:], fit_result, 
                               save_dir+'/images/fwhm_fit/', str(i).zfill(4)+'fit_im0.png')
    
    return result

def makeGaussian(amplitude,x0=256,y0=512,sigx=4,sigy=4,CCDsize=[1024,1024],b=0):
    x = np.arange(0,CCDsize[0],1)
    y = np.arange(0,CCDsize[1],1)
    #print(x,y)
    x, y = np.meshgrid(x,y)
    #print(xy)
    GX = (x-x0)**2 / (2*sigx**2)
    GY = (y-y0)**2 / (2*sigy**2)
    Gaussian2D = amplitude * np.exp(-(GX+GY)) + b
    return Gaussian2D

def psf_night(save_dir,amplitude,night,Xc=256,Yc=512,sigma_fwhm=4.,CCDsize=[1024,1024],
              rot=0,b=0,save_plot=True,show=False,show_time=True):
    '''
    Create the simulated psf for a night based on real data for the FWHM
    
    CCDsize = [1024,1024]
    x0,y0,a,sigx,sigy,rot,b = 256,512,8000,4,4,0,0
    '''
#     print(Xc,Yc)
    start = Time.time()
#     sigma_fwhm = fwhmfit.fwhm_xo2b(len(amplitude),save_dir+'fwhm_data.csv')
    #use 
    sigx, sigy = np.ones(len(amplitude))*sigma_fwhm,np.ones(len(amplitude))*sigma_fwhm
#     print(len(sigx),len(sigy))
    
    #creating and alocate space for some variables
    rawflux,signal_ = np.zeros(size), np.zeros(size)
            
    df_rawflux,df_signal = pd.DataFrame(), pd.DataFrame()
    ica_norm_df, img_lines, img_columns = pd.DataFrame(), np.zeros(size), np.zeros(size)
    midpointX, midpointY = np.zeros(size),np.zeros(size)
    amplitude_fit = np.zeros(size)
            
    #width of the Gaussian PSF model
    wid = np.zeros(size)
    #fwhm retrived from the model
    fwhm = np.zeros(size)
    fwhm_norm = np.zeros(size)
    
    #creating the tree directories to save files
    if not os.path.exists(save_dir+night):
        os.makedirs(save_dir+night)
    if not os.path.exists(save_dir+night+'/snapshot'):
        os.makedirs(save_dir+night+'/snapshot')
    if not os.path.exists(save_dir+night+'/images/normPSF/'):
        os.makedirs(save_dir+night+'/images/normPSF/')
    if not os.path.exists(save_dir+night+'/images/fwhm_fit/'):
        os.makedirs(save_dir+night+'/images/fwhm_fit/')
    
    #creating observation
    for i in range(len(amplitude)):
        xx,yy = drift_model(save_dir_orig,Xc,Yc,len(amplitude))
        
        img = makeGaussian(amplitude=amplitude[i])
#         print(img)
#         break
#         img, star, pars_psf, area = fwhmfit.exposition(int(xx[i]),int(yy[i]),amplitude[i],int(sigx[i]),int(sigy[i]),rot,b,CCDsize,show = False)
#         img = img.data
        
        if save_plot == True:
            plt.figure(figsize=(11,9))
            plt.imshow(img,cmap=plt.cm.coolwarm)
            plt.xlim(x0-10,x0+10) #center x0 = 256
            plt.ylim(y0-10,y0+10) #center y0 = 512
            plt.colorbar(label='# counts')
            plt.xlabel('x [px]')
            plt.ylabel('y [px]')
            plt.xticks(np.linspace(x0-10,x0+10,5))
            plt.yticks(np.linspace(y0-10,y0+10,5))
            plt.savefig(save_dir+night+'/snapshot/'+str(i).zfill(5)+'.png')
            if show == True:
                plt.show()
                plt.close()
            else:
                plt.close()

        #FItting the center of the PSF
        delta = 14
        img = img[int(512-delta):int(512+delta),int(256-delta):int(256+delta)]
        centerXY = fwhmfit.centerPSF(img)
        fit_result =  fwhmfit.psf_model(np.arange(0,img.shape[0],1),
                                        img[int(centerXY[0]),:],model='gaussian', 
                                        amp=img[int(centerXY[0]),:][img[int(centerXY[0]),:].argmax()]/2.,
                                        cen = centerXY[1], wid = 1.,show=False)
        if save_plot == True:
            fwhmfit.psf_model_plot(np.arange(0,img.shape[0],1),img[int(centerXY[0]),:], fit_result, 
                               save_dir+night+'/images/fwhm_fit/', str(i).zfill(4)+'fit_im0.png')
        
        # Apperture Photometry using Astropy
        radius = 2.5* np.sqrt(np.sqrt(fit_result.params['wid']/2.))
#         print(radius)
        aperture = CircularAperture(centerXY,r=2)
        flux = aperture_photometry(img,aperture)
        rawflux[i] = flux['aperture_sum']
#         break
        
        #Saving the results
        midpointX[i] = centerXY[1]
        midpointY[i] = fit_result.params['cen']
        amplitude_fit[i] = fit_result.params['amp']
        wid[i] = fit_result.params['wid']
        fwhm[i] = 2.355 * np.sqrt(fit_result.params['wid']/2.)
#         rawflux[i] = fwhmfit.flux_aperture(img,centerXY[0],centerXY[1],7)
#         break
    
    if show_time == True:
        print('Total time = ',round(abs(Time.time()-start)/60.,2))
        
    #creating dicctionary for results
    result ={
        'rawflux':rawflux,
        'midpointX':midpointX,
        'midpointY':midpointY,
        'amplitude_fit':amplitude_fit,
        'wid':wid,
        'fwhm':fwhm,
        'centerXY':centerXY
    }
    return result

