import os
import glob
import numpy as np
import astropy.io.fits as fits
# from ExoSetupTaskParameters import * #PyIraf configuration
# import pyraf as iraf
# from loginIraf import * #loading login.cl parameters for iraf
from pandas import DataFrame
from scipy import misc #work with autocorrelation in images

def imcombine(images,result):
    '''
    Combine the images using the median value as the default.
    images: list of string with the names of the images to be combined.
    result: string, name of the combined image
    '''
    im_array = []
    for i in range(len(images)):
        im = fits.getdata(images[i])
        im_array.append(np.array(im, dtype='Float64'))
    image_array = np.median(im_array, axis=0)
    hdu_images = fits.PrimaryHDU(image_array)
    hdulist_images = fits.HDUList([hdu_images])
    hdulist_images.writeto(result)

def imarith(spc_image,images,procedure ='subtracted',sulfix='b'):
    '''
    Divided or subtracted one specfic image from a group of images.
    spc_image: array-like, matrix to be subtracted or divided from the group of images
    images: list of strings, name of the images to subtract the matrix spc_image
    '''
    if procedure == 'subtracted':
        print('\n Subtracting spc_image from the images and creating sulfix*images.... \n')
        for i in range(len(images)):
            im = fits.getdata(images[i])
            im = np.array(im,dtype='Float64')
            bimage = im - spc_image
            hdu_images = fits.PrimaryHDU(bimage)
            hdulist_images = fits.HDUList([hdu_images])
            hdulist_images.writeto(sulfix+images[i])
    if procedure == 'divided':
        print('\n Dividing the spc_image from the images and creating sulfix*images.... \n')
        for i in range(len(images)):
            im = fits.getdata(images[i])
            im = np.array(im,dtype='Float64')
            bimage = im / spc_image
            hdu_images = fits.PrimaryHDU(bimage)
            hdulist_images = fits.HDUList([hdu_images])
            hdulist_images.writeto(sulfix+images[i])

def masterbias(bias_sulfix,data_path,save_path):
    """
    Obtain the masterbias.fits image.
    ___
    Input:
    For obtain this parameters, use the input_info function.

    data_path: string, path where are the images data.
    save_path: string, path where will save all reduced images.
    input_file: dict, with information describe in the YAML file.

    Output:
    It is possible that the function return some of these values:

    0. Create the masterbias image on the save_path.
    1. It do not create the masterbias image, because of some error
    ___
    """
    #Set original directory
    original_path = os.getcwd()
    #Change your directory to data diretory
    os.chdir(data_path)
    #list all bias images
    bias = glob.glob(bias_sulfix+'*.fits')
    print('Loading bias images \nTotal of bias files = ',len(bias),'\nFiles = \n')
    print(bias)
    print('\nCreating superbias \n')
    #if save_path exist, continue; if not, create.
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    #copy bias images to save_path
    os.system('cp '+bias_sulfix+'*.fits '+save_path)
    #change to sabe_path
    os.chdir(save_path)
    #verify if previous superbias exist
    if os.path.isfile('superbias.fits') == True:
        os.system('rm superbias.fits')
    # --------------------------------------------------------------------------
    # --- Using only with a few bias images
    #create the list of bias images
    #bias_list = string.join(bias,',')
    #combine the bias image and create the superbias
    #iraf.imcombine(bias_list,'superbias.fits')
    #iraf.imstat('superbias.fits')
    # --------------------------------------------------------------------------

    #Using numpy package to take the mean value of bias images
    #Problem: does not include the superbias header in this version
    # bias_array = []
    # for i in range(len(bias)):
    #     image = fits.getdata(bias[i])
    #     bias_array.append(np.array(image,dtype='Float64'))
    # superbias_array = np.median(bias_array,axis=0)
    # hdu_superbias = fits.PrimaryHDU(superbias_array)
    # hdulist_superbias = fits.HDUList([hdu_superbias])
    # hdulist_superbias.writeto('superbias.fits')
    imcombine(bias,'superbias.fits')
    #clean previos bias files
    print('\n Cleaning bias*.fits images ....\n')
    os.system('rm '+bias_sulfix+'.fits')
    print('\n.... done.')
    #print output
    #test of outpu value
    #os.remove('superbias.fits')
    #Verify if the image was created:
    output = glob.glob('superbias*.fits')
    if len(output) != 0:
        output = 0
    else:
        output = 1
    #Return to original directory
    os.chdir(original_path)
    #END of the masterbias reduction messsage
    print('\nsuperbias.fits created!\n')
    print('\nEND of superbias reduction!\n')
    #obtain the value of return
    if output == 1:
        print('!!! ERROR/WARNING !!!')
        print('Check if the superbias was created or if there is more than one superbias image.')
    return output

def masterflat(flat_sulfix,data_path,save_path):
    """
    Obtain the masterflat image for calibration.
    ___
    INPUT:
    For obtain this parameters, use the input_info function.

    data_path: string, path where are the images data.
    save_path: string, path where will save all reduced images.
    input_file: dict, with information describe in the YAML file.

    OUTPUT:
    It is possible that the function return some of these values:

    0. Create the masterflat image on the save_path.
    1. It do not create the masterflat image, because of some erros.
    """
    #set original directory
    original_path = os.getcwd()
    #Change your directory to data diretory
    os.chdir(data_path)
    #list all flat images
    flat = glob.glob(flat_sulfix+'*.fits')
    # Check if exist flat images
    if not flat:
        print('there is no flat images to be reduced')
        os.chdir(original_path)
        return 1

    print('Loading flat images \nTotal of flat files = ', len(flat), '\nFiles = \n')
    print(flat)
    #if save_path exist, continue; if not, create.
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    #create a list of bias images and copy images to save_path
    os.system('cp '+flat_sulfix+'*.fits '+save_path)
    #creating the names of flat with bias subctracted
    bflat = []
    for i in flat:
        bflat.append('B'+i)
    print('\n Names os flat images with bias subtracted: \n \n', bflat)
    #change for save_path directory
    os.chdir(save_path)
    #verify if previous superbias exist
    if os.path.isfile('superflat.fits') == True:
        os.system('rm superflat.fits')
    #verify if exits previous bflat*.fits files and remove then.
    for i in bflat:
        if os.path.isfile(i) == True:
            os.system('rm -f '+i)
    print('\nCreating superflat .... \n')
    #************************************************************************************************************
    #******************  BEGIN Pyraf reduction ( IRAF reduction in Python) **************************************
    #************************************************************************************************************
    #create the list of flat images  and bflat images
    #flat = string.join(flat,',')
    #bflat = string.join(bflat,',')
    # print '\n Subtracting bias from flat images and creating bflat images.... \n'
    # #iraf.imarith()
    # for i in range(len(flat)):
    #     iraf.imarith(flat[i],'-','superbias.fits',bflat[i])
    #     #print statistics from bflat*.fits images
    #     iraf.imstat(bflat[i])
    # print '\n .... done \n'
    #************************************************************************************************************
    #********************  END Pyraf reduction ( IRAF reduction in Python) **************************************
    #************************************************************************************************************

    #************************************************************************************************************
    #******************  BEGIN Python Numpy Redction with Astropy.io.fits  **************************************
    #************************************************************************************************************
    masterbias_spc = fits.getdata('superbias.fits',header=False)
    masterbias_spc = np.array(masterbias_spc,dtype='Float64')
    imarith(masterbias_spc,flat,procedure ='subtracted',sulfix='b')
    #************************************************************************************************************
    #********************  END Python Numpy Redction with Astropy.io.fits  **************************************
    #************************************************************************************************************

    #clean previos flat*.fits files
    print('\n Clean flat*.fits images .... \n')
    os.system('rm '+flat_sulfix+'*.fits')
    print('\n .... done. \n')
    #normalizing each flat
    print('\nNormalizing each flat ....\n')
    #checking if mean from numpy is the same from your bflat images using imstat
    #take the mean of each bflat image
    bflat_mean = np.zeros(len(bflat))
    for i in range(len(bflat)):
        image = fits.getdata(bflat[i])
        image = np.array(image, dtype='Float64')
        bflat_mean[i] = round(np.mean(image))
    image = 0 #clean image allocate to this variable
    print('The mean of each bflat image, respectivaly ...')
    print(bflat_mean)
    #creating the names of bflat images after the normalization:
    abflat = []
    for i in bflat:
        abflat.append('A'+i)
    print('\n Names os bflat images with bias subtracted and normalizad: \n \n', abflat)
    #verify if exist previous ABflat*.fits images and remove then.
    for i in abflat:
        if os.path.isfile(i) == True:
            os.system('rm -f '+i)
    
    #************************************************************************************************************
    #******************  BEGIN Pyraf reduction ( IRAF reduction in Python) **************************************
    #************************************************************************************************************
    # for i in range(len(abflat)):
    #     iraf.imarith(bflat[i],'/',bflat_mean[i],abflat[i])
    # print '\n.... done!\n'
    # # print '\n Cleaning bflat*.fits images ....\n'
    # # os.system('rm Bflat*.fits')
    # print '\n.... done.\n'
    # print 'Statistics of the abflat*.fits images .... \n'
    # for i in range(len(abflat)):
    #     iraf.imstat(abflat[i])
    # print '\n Combining abflat images ....\n'
    # #***
    # ablist = string.join(abflat,',')
    # iraf.imcombine(ablist,'superflat.fits')
    #************************************************************************************************************
    #********************  END Pyraf reduction ( IRAF reduction in Python) **************************************
    #************************************************************************************************************

    #************************************************************************************************************
    #******************  BEGIN Python Numpy Redction with Astropy.io.fits  **************************************
    #************************************************************************************************************
    for i in range(len(abflat)):
        print('working on '+abflat[i])
        bflat_spc = fits.getdata(bflat[i],header=False)
        bflat_spc = np.array(bflat_spc,dtype='Float64')
        bflat_spc = np.ones(bflat_spc.shape)*bflat_mean[i]
        print('Matrix of normalization = ',bflat_spc)
        imarith(bflat_spc,[bflat[i]],procedure ='divided',sulfix='A')
        abflat_im = fits.getdata(abflat[i],header=False)
        abflat_im = np.array(abflat_im,dtype='Float64')
        print('Mean',round(np.mean(abflat_im)),' STD = ',round(np.std(abflat_im)))
        print('\n.... done!\n')
    #************************************************************************************************************
    #********************  END Python Numpy Redction with Astropy.io.fits  **************************************
    #************************************************************************************************************
    print('\n Combining abflat images ....\n')
    #change how import flat files
    #usning the abflat list of flat files We will create a pandas python dataframe
    # ablist = DataFrame(abflat)
    # ablist.columns=['flat_files']
    # ablist.to_csv('flat_list',index_label=False,index=False,header=False)

    #************************************************************************************************************
    #******************  BEGIN Pyraf reduction ( IRAF reduction in Python) **************************************
    #************************************************************************************************************
    # #combine all flat images
    # iraf.imcombine('@flat_list','superflat.fits')
    # iraf.imstat('superflat.fits')
    # print '\n .... done. \n'
    # #************************************************************************************************************
    #********************  END Pyraf reduction ( IRAF reduction in Python) **************************************
    #************************************************************************************************************

    #************************************************************************************************************
    #******************  BEGIN Python Numpy Redction with Astropy.io.fits  **************************************
    #************************************************************************************************************
    imcombine(abflat,'superflat.fits')
    #************************************************************************************************************
    #********************  END Python Numpy Redction with Astropy.io.fits  **************************************
    #************************************************************************************************************

    # print '\nCleaning ABflat*.fits images ....\n'
    # os.system('rm ABflat*.fits')
    print('\n.... done!')
    #Verify if the image was created:
    output = glob.glob('superflat*.fits')
    if len(output) != 0:
        output = 0
    else:
        output = 1
    #Return to original directory
    os.chdir(original_path)
    #last mensage
    print('\n MASTERFLAT.FITS created! \n')
    print('\n END of Data Reduction for create a masterflat.fits file. \n')
    #obtain the value of return
    if output == 1:
        print('!!! ERROR/WARNING !!!')
        print('Check if the superbias was created or if there is more than one superbias image.')
    return output


def science_images(data_path,save_path,image_sulfix,abimages=None):
    '''
    Calibrate science images with masterflat (or superflat) and masterbias (or superbias) images.
    '''
    #set original directory
    original_path = os.getcwd()
    os.chdir(data_path)
    sci_images = np.sort(glob.glob(image_sulfix+'*.fits'))
    print('science images = \n',sci_images)
    #if save_path exist, continue; if not, create.
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    #create a list of bias images and copy images to save_path
    print('\nCopy science images to save_path directory to main reduction: ....')
    os.system('cp '+image_sulfix+'*.fits '+save_path)
    print('\n .... done. \n')
    #change to save_path
    os.chdir(save_path)
    ######################################################################################################
    #################### Creating the names of the images  ###############################################
    ######################################################################################################
    b_sci_images,ab_sci_images = [],[]
    for i in sci_images:
        ab_sci_images.append('AB'+i)
        b_sci_images.append('B'+i)
    ######################################################################################################
    #################### Removing the bias image from the science images #################################
    ######################################################################################################
    #create the names for exoplanet science mages with bias subtracted
    for i in sci_images:
        #verify if previous superbias exist
        if os.path.isfile('B'+i) == True:
            os.system('rm B'+i)
    bias_im = fits.getdata('superbias.fits',header=False)
    bias_im = np.array(bias_im,dtype='Float64')
    imarith(bias_im,sci_images,procedure ='subtracted',sulfix='B')
    os.system('rm '+image_sulfix+'*.fits')

    ######################################################################################################
    #################### Removing the flat image from the science images #################################
    ######################################################################################################
    for i in sci_images:
        #verify if previous superbias exist
        if os.path.isfile('AB'+i) == True:
            os.system('rm AB'+i)
    flat_im = fits.getdata('superflat.fits',header=False)
    flat_im = np.array(flat_im,dtype='Float64')
    imarith(flat_im,b_sci_images,procedure ='divided',sulfix='A')
    ######################################################################################################
    ######################################### science images #############################################
    ######################################################################################################
    #clean up the folder from previous images
    os.system('rm B'+image_sulfix+'*.fits')
    print('Science images corrected by bias and flatfield: \n',ab_sci_images)
    os.chdir(original_path)
    return 0

def phase_correlation(a, b):
    '''
    Autocorrelation between two images, a and b, pixel a pixel.
    '''
    G_a = np.fft.fft2(a)
    G_b = np.fft.fft2(b)
    conj_b = np.ma.conjugate(G_b)
    R = G_a*conj_b
    R /= np.absolute(R)
    r = np.fft.ifft2(R).real
    return r