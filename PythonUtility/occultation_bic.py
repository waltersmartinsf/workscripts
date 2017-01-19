import numpy as np
from scipy.stats import chisquare
#You can download at: http://astro.uchicago.edu/~kreidberg/
from occultquad import occultquad
#mpfit function in python from the package mpyfit
#You can download at: https://github.com/waltersmartinsf/mpyfit
import mpyfit
from pandas import read_csv, DataFrame,concat
import glob
import time
import sys
import os

"""
BIC Procedure for Primary Transit of Exoplanets
"""
#BAR Progress function to visualize the progress status:
def update_progress(progress):
    """
    Progress Bar to visualize the status of a procedure
    ___
    INPUT:
    progress: percent of the data

    ___
    Example:
    print ""
    print "progress : 0->1"
    for i in range(100):
        time.sleep(0.1)
        update_progress(i/100.0)
    """
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

def sysparams(path):
    """
    Obtain the sysparam*.txt information of the exoplanet.

    INPUT

    path: the sysparams string path in the system:

    Example:

    a_Rs, Rp_Rs, radstar, i, a, e, omega, epoch, period, limbB1, limbB2 = sysparams('/home/user/exoplanet_a/filter_B/sysparams_B.txt')
    """
    a_Rs, Rp_Rs, radstar, i, a, e, omega, epoch, period, limbB1, limbB2 = np.loadtxt(path)
    return a_Rs, Rp_Rs, radstar, i, a, e, omega, epoch, period, limbB1, limbB2

def read_files(workpath,hjd_path,airmass_path):
    """
    workpath: string, path to allflux* and error_allflux* files
    ___
    Example:
    #hjd is the ascii file with Julian Dates
    hjd_path = '/home/user/planet/Results/hjd'
    ___

    OUTPUT
    rawflux: relative flux, normalized
    eflux: associated error flux
    """
    #save the current directory
    original_path = os.getcwd()
    #change for the work diretory
    os.chdir(workpath)

    #Reading data
    # hjd = read_csv(hjd_path,names=['hjd']) #julian date for each point in the observation night
    # hjd = hjd.hjd.values #our JD dates from ExoDPRL
    hjd = np.loadtxt(hjd_path)
    x,y,am = np.loadtxt(airmass_path,unpack=True) #airmass
    files = np.sort(glob.glob('allflux*')) #create the variable allflux with the names of our allflux files in the directory
    erro_files = np.sort(glob.glob('error_allflux*'))
    files = DataFrame([files,erro_files]).T
    files.columns = ['files','Err_files']
    filetot = len(files)
    print 'Total of files = ',filetot
    #Verify the standart deviation of our data
    scatter = np.zeros(filetot)
    for i in range(filetot):
        #print('Working @'+files.files[i])
        #print('... Read '+files.files[i]+' ...')
        exoplanet_data = read_csv(files.files[i],sep=' ',header=None)#,names=['xo2n','xo2s'])
        #print('... done.')
        #standart deviation in n data from allflux
        scatter[i] = exoplanet_data[0].std()
        update_progress((i+1.)/filetot)
    files[2] = scatter #include the scatter datain the files dataframe structure
    files.columns = ['files','Err_files','scatter']
    id_min = files['scatter'].argmin() #index of the min scatter file
    id_max = files['scatter'].argmax() #index for the maximum scatter file
    print 'The smallest scatter is: '+str(files.scatter[id_min])
    print 'Which is file: '+files.files[id_min]
    #Obtain the rawflux

    print('Working @'+files.files[id_min]+' that is the min scatter')
    print('... Read '+files.files[id_min]+' ...')
    exoplanet_data = read_csv(files.files[id_min],sep=' ',header=None)#,names=['xo2n','xo2s'])
    print('... done.')
    print('... Read '+files.Err_files[id_min]+' ...')
    Err_exoplanet_data = read_csv(files.Err_files[id_min],sep=' ',header=None)#,names=['Err_xo2n','Err_xo2s'])
    print('... done')

    #obtain the comparitive flux:
    rows_number,columns_number = exoplanet_data.shape
    #print exoplanet_data.shape
    
    #obtain the comparitive flux:
    flux = DataFrame()
    eflux = DataFrame()
    for j in range(int(columns_number)):
        if j != 0:
            _flux = exoplanet_data[0].values/exoplanet_data[j].values
            part1 = (Err_exoplanet_data[0]/exoplanet_data[0])**2
            part2 = (Err_exoplanet_data[j]/exoplanet_data[j])**2
            _eflux = _flux*np.sqrt(np.array(part1.values + part2.values))
            _flux = DataFrame(_flux,columns=[str(j-1)])
            _eflux = DataFrame(_eflux,columns=[str(j-1)])
            frames = [flux,_flux]
            flux = concat(frames,axis=1)
            frames = [eflux,_eflux]
            eflux = concat(frames,axis=1)
    rawflux = flux
    count_flux = exoplanet_data[0].values
    count_eflux = Err_exoplanet_data[0].values
    os.chdir(original_path)
    return hjd, rawflux, eflux, am, count_flux,count_eflux, exoplanet_data, Err_exoplanet_data

#*******************************************************************************************************************************
#*******************************************************************************************************************************
# LIGTHCURVE MODEL

def delta(phase,inc, ecc = 0, omega=0):
    """
    Compute the distance center-to-center between planet and host star.
    ___

    INPUT:

    phase: orbital phase in radian
    inc: inclination of the system in radian

    OPTIONAL INPUT:

    ecc:
    omega:

    //
    OUTPUT:

    distance center-to-center, double-float number.
    ___


    """
    phase = 2*np.pi*phase
    if ecc == 0 and omega == 0:
        delta = np.sqrt(1-(np.cos(phase)**2)*(np.sin(inc)**2))
    else:
        delta = (1.-ecc**2.)/(1.-ecc*np.sin(phase-omega))* np.sqrt((1.-(np.cos(phase))**2.*(np.sin(inc))**2))

    return delta

def limb_darkening(limbb1, limbb2, mu):
    """
    Compute the limb-darkening at the edge of the star.
    This function calcilate the lim darkening with a quadratic function,
    following the ideas of Mandel & Agol (2002), and Winn (2010).
    ___
    mu = cos(incident angle)
    limb1, limbb2 host-star limb darkening parameters (from sysparams*.txt)
    """
    icalc = 1.
    icalc = icalc - limbb1*(1-mu) - limbb2*((1-mu)**2.)
    return icalc

def occultation_fast(z,w,show=False):
    """
    Obtain the lightcurve without limb darkening.
    ___
    INPUT:
    z: impact parameter in units of rs (list)
    w: occulting star size in units of rs (list)

    OUTPUT:
    muo1 : fraction of flux at each b0 for a uniform source
    """
    if isinstance(z,(int,float)) == True: #change int or float to np.ndarray
        z = [z]
    if isinstance(w,(int,float)) == True:
        w = [w]

    z = np.array(z)
    w = np.array(w)

    if show == True:
        print('Impact parameters, z: ',z)

    if abs(w.all() - 0.5) < 1.E-3: #line to garantee that w is evaluate for all possibilities os the array.
        w = 0.5
    if show == True:
        print('Occultating size star, w: ',w)

    nz = len(z)
    muo1 = np.zeros(nz) #create the return array of the problem
    indx = np.where(z > 1.+ w ) #first condition: impact parameter greater than occulting star size
    #If this condition is true, then where this condition is true, the fraction of flux is equal to 1.
    indx = indx[0] #choose the first index --> list of index that the condition of where is TRUE
    if indx.size == 0:
        if show == True:
            print('(1) Any argument inside the impact parameters, z, satisfy the condition: z > 1.+w')
    else:
        if show == True:
            print('(1) The impact parameters satisfy the condition: z <= 1.+w')
        muo1[indx] = 1.
    ####
    #Second condition:
    #impact parameter z >=  |1. - occultation star size|  AND impact parameter z <= (1. + occultation star size)
    indx = np.where((z >= abs(1.-w)) & (z <= 1.+ w))
    indx = indx[0]
    if indx.size == 0:
        if show == True:
            print('(2) Any argument inside the impact parameters, z, satisfy the condition: z >= |1.-w| and z<= 1.+w')
    else:
        if show == True:
            print('(2) Impact parameters, z, satify the condition: z < |1.- w| and z > 1. + w')
        zt = z[indx]
        xt=(1.-w**2+zt**2.)/2./zt

        kap1 = np.zeros(len(xt))
        for i in range(len(xt)):
            if xt[i] < 1. :
                kap1[i] = np.arccos(xt[i])
            #xt[i] >= 1. :
            else:
                kap1[i] = np.arccos(1.)

        xt=(w**2+zt**2-1.)/2./w/zt

        kap0 = np.zeros(len(xt))
        for i in range(len(xt)):
            if xt[i] < 1.:
                kap0[i] = np.arccos(xt[i])
            else: # xt[i]>= 1.
                kap0[i] = np.arccos(1.)

        lambdae=w**2*kap0+kap1
        xt=4.*zt**2-(1.+zt**2-w**2)**2

        for i in range(len(lambdae)):
            if xt[i] >= 0:
                lambdae[i] = (lambdae[i] - 0.5 * np.sqrt(xt[i]))/np.pi
            else:
                lambdae[i] = lambdae[i]/np.pi
        muo1[indx]=1.-lambdae
    ####
    ### Third condition
    indx = np.where(z < 1.-w)
    indx = indx[0]
    if indx.size == 0:
        if show == True:
            print('(3) Any argument inside the impact parameters, z, satisfy the condition: z < 1.- w')
    else:
        if show == True:
            print('(3) Impact parameters,z, satisfy the condition: z >= 1. - w')
        muo1[indx]=1.-w**2
    return muo1

#*******************************************************************************
#*******************************************************************************
def simple_model(JD, period, inc, a_Rs, limbB1, limbB2,startpar0, startpar1):
    phase = (JD - startpar1)/period
    distance_vector = delta(phase,inc) * a_Rs
    #model = occultation_fn(distance_vector,startpar0,limbB1,limbB2,show=False)
    model = occultquad(distance_vector, limbB1, limbB2, startpar0, len(JD))

    return model

def model_am_exp(hjd, period, inc, a_Rs, limbB1, limbB2,startpar0,startpar1,startpar2,startpar3):
    model = simple_model(hjd, period, inc, a_Rs, limbB1, limbB2,startpar0, startpar1)
    model_am = model * startpar2 * np.exp(-1. * startpar3 * am) #multiply light curve by factor x exponential
    return model_am

def residuals_am_exp(params,args): #residuals function for mpfit
    RpRs = params[0]
    Tc = params[1]
    mu1 = params[2]
    mu2 = params[3]
    hjd, data, eps_data, period, inc, a_Rs, limbB1, limbB2 = args
    model = model_am_exp(hjd, period, inc, a_Rs, limbB1, limbB2,RpRs,Tc,mu1,mu2)
    return (data-model)/eps_data

def model_am_linear(hjd, period, inc, a_Rs, limbB1, limbB2,startpar0,startpar1,startpar2,startpar3):
    model = simple_model(hjd, period, inc, a_Rs, limbB1, limbB2,startpar0, startpar1)
    model_am = model * (startpar2 * am + startpar3) #multiply light curve by factor x exponential
    return model_am

def residuals_linear(params,args): #residuals function for mpfit
    RpRs = params[0]
    Tc = params[1]
    mu1 = params[2]
    mu2 = params[3]
    hjd, data, eps_data, period, inc, a_Rs, limbB1, limbB2 = args
    model = model_am_linear(hjd, period, inc, a_Rs, limbB1, limbB2,RpRs,Tc,mu1,mu2)
    return (data-model)/eps_data

def model_am_2deg(hjd, period, inc, a_Rs, limbB1, limbB2,startpar0,startpar1,startpar2,startpar3,startpar4):
    model = simple_model(hjd, period, inc, a_Rs, limbB1, limbB2,startpar0, startpar1)
    model_am = model * (startpar2 + startpar3*am + startpar4*am) #multiply light curve by factor x exponential
    return model_am

def residuals_2deg_mpfit(params,args): #residuals function for mpfit
    RpRs = params[0]
    Tc = params[1]
    mu1 = params[2]
    mu2 = params[3]
    mu3 = params[4]
    hjd, data, eps_data, period, inc, a_Rs, limbB1, limbB2 = args
    model = model_am_2deg(hjd, period, inc, a_Rs, limbB1, limbB2,RpRs,Tc,mu1,mu2,mu3)
    return (data-model)/eps_data

def BIC(nfree,bestnorm, nflux):
    """
    Obtain the BIC value for the model.
    ___
    INPUT:
    nfre : number of free parameters
    bestnorm:
    nflux: number of points in the data


    OUTPUT:
    BIC: BIC number
    """
    if isinstance(nfree,(float,int)) == True:
        bic = nfree * np.log(nflux) + bestnorm
    else:
        bic = np.zeros(len(nfree))
        for i in range(len(nfree)):
            bic[i] = nfree[i] * np.log(nflux) + bestnorm[i]
    return bic
# 
# def model_airmassfit(hjd, am, rawflux, eflux, limbB1, limbB2, inc, period, a_Rs, Rp_Rs, show=False):
#     """
#     Return the bestfit model for the lightcurve using 4 models of airmass correction:
#     1. model with no airmass correction
#     2. model with exponential airmass correction
#     3. model with linear airmass correction
#     4, model with 2deg polynomial airmass correction
#     ___
#     INPUT:
#
#     hjd:
#     am:
#     rawflux:
#     eflux:
#     limbB1:
#     limbB2:
#     inc:
#     period:
#     a_Rs:
#     startpar:
#
#     OUTPUT:
#     result: dataframe structure with besfit values for each model, the errors and BIC values.
#     phase: from the bestfit model
#     lc: lightcurve from the bestfit model
#     """
#     # Model 1: no airmass correction
#     startpar = [Rp_Rs, np.mean(hjd), 1., 0.]
#     PARINFO = [{'value':Rp_Rs,'limits':(0,1.)}, {'value':np.mean(hjd)}, {'value':1.}, {'value':0.,'fixed':True}]
#     pfit1, results1 = mpyfit.fit(residuals_am_exp, startpar, args = (hjd,rawflux,eflux, period, inc, a_Rs, limbB1, limbB2), parinfo=PARINFO)
#     model1 = model_am_exp(hjd, period, inc, a_Rs, limbB1, limbB2,pfit1[0],pfit1[1],pfit1[2],pfit1[3])
#     phase1 = (hjd - pfit1[1])/period
#     if show == True:
#         print '...'
#         print 'Model 1: no airmass correction'
#         print 'bestfit values = ',pfit1
#         print 'error = ', results1['parerrors']
#         print 'bestnorm1 = ', results1['bestnorm']
#         print 'chi-square, scipy routine = ',chisquare(rawflux, model1)
#     #Model 2: exponential airmass correction
#     PARINFO = [{'value':Rp_Rs,'limits':(0,1.)}, {'value':np.mean(hjd)}, {'value':1.}, {'value':0.,'fixed':False}]
#     pfit2, results2 = mpyfit.fit(residuals_am_exp, startpar, args = (hjd,rawflux,eflux, period, inc, a_Rs, limbB1, limbB2), parinfo=PARINFO)
#     model2 = model_am_exp(hjd, period, inc, a_Rs, limbB1, limbB2,pfit2[0],pfit2[1],pfit2[2],pfit2[3])
#     phase2 = (hjd - pfit2[1])/period
#     if show == True:
#         print '...'
#         print 'Model 2: exponential airmass correction'
#         print 'bestfit values = ',pfit2
#         print 'error = ', results2['parerrors']
#         print 'bestnorm1 = ', results2['bestnorm']
#         print 'chi-square, scipy routine = ',chisquare(rawflux, model2)
#     #Model 3: linear airmass correction
#     PARINFO = [{'value':Rp_Rs,'limits':(0,1.)},{'value':np.mean(hjd)},{'value':1.}, {'value':0.,'fixed':False}]
#     pfit3, results3 = mpyfit.fit(residuals_linear, startpar, args = (hjd,rawflux,eflux, period, inc, a_Rs, limbB1, limbB2), parinfo=PARINFO)
#     model3 = model_am_linear(hjd, period, inc, a_Rs, limbB1, limbB2,pfit3[0],pfit3[1],pfit3[2],pfit3[3])
#     phase3 = (hjd - pfit3[1])/period
#     if show == True:
#         print '...'
#         print 'Model 3: linear airmass correction'
#         print 'bestfit values = ',pfit3
#         print 'error = ', results3['parerrors']
#         print 'bestnorm1 = ', results3['bestnorm']
#         print 'chi-square, scipy routine = ',chisquare(rawflux, model3)
#     #Model 4: 2deg polynomial airmss correction
#     PARINFO = [{'value':Rp_Rs,'limits':(0,1.)},{'value':np.mean(hjd)},{'value':1.},{'value':0.},{'value':0.}]
#     pstart = [Rp_Rs,np.mean(hjd),1.,0.,0.]
#     pfit4, results4 = mpyfit.fit(residuals_2deg_mpfit, pstart, args = (hjd,rawflux,eflux, period, inc, a_Rs, limbB1, limbB2), parinfo=PARINFO)
#     model4 = model_am_2deg(hjd, period, inc, a_Rs, limbB1, limbB2,pfit4[0],pfit4[1],pfit4[2],pfit4[3],pfit4[4])
#     phase4 = (hjd - pfit4[1])/period
#     if show == True:
#         print '...'
#         print 'Model 4: 2deg poly airmass correction'
#         print 'bestfit values = ',pfit4
#         print 'error = ', results4['parerrors']
#         print 'bestnorm1 = ', results4['bestnorm']
#         print 'chi-square, scipy routine = ',chisquare(rawflux, model4)
#     #Obtain BIC values:
#     #Let's create our fit file and our best BIC
#     BICarray = ['none', 'exponential', 'linear','2nd_deg_poly']
#     nfree = [3,4,4,5]
#     bestnorm = [results1['bestnorm'],results2['bestnorm'],results3['bestnorm'],results4['bestnorm']]
#     bic = BIC(nfree,bestnorm,len(rawflux))
#     RpRs = [pfit1[0], pfit2[0], pfit3[0], pfit4[0]]
#     Tc   = [pfit1[1], pfit2[1], pfit3[1], pfit4[1]]
#     a    = [pfit1[2], pfit2[2], pfit3[2], pfit4[2]]
#     b    = [pfit1[3], pfit2[3], pfit3[3], pfit4[3]]
#     c    = ['Nan','Nan','Nan',pfit4[4]]
#     error1 = [results1['parerrors'][0], results2['parerrors'][0], results3['parerrors'][0], results4['parerrors'][0]]
#     error2 = [results1['parerrors'][1], results2['parerrors'][1], results3['parerrors'][1], results4['parerrors'][1]]
#     error3 = [results1['parerrors'][2], results2['parerrors'][2], results3['parerrors'][2], results4['parerrors'][2]]
#     error4 = [results1['parerrors'][3], results2['parerrors'][3], results3['parerrors'][3], results4['parerrors'][3]]
#     error5 = ['Nan','Nan','Nan', results4['parerrors'][0]]
#     result = DataFrame([BICarray,list(bic),RpRs,error1,Tc,error2,a,error3,b,error4,c,error5]).T
#     result.columns=['Model','BIC','RpRs','eRpRs','Tc','eTc','a','ea','b','eb','c','ec']
#     if show == True:
#         print '... Results:'
#         print result
#         print 'The best model is: ',result.Model[result.BIC == result.BIC.min()]
#         print 'with the BIC = ',result.BIC.min()
#     #Saving the bestfit transit image:
#     bestfit = np.where(result.BIC == result.BIC.min())
#     indx = bestfit[0][0]
#     if indx == 0:
#         lc = model1
#         phase = phase1
#     if indx == 1:
#         lc = model2
#         phase = phase2
#     if indx == 2:
#         lc = model3
#         phase = phase3
#     if indx == 3:
#         lc = model4
#         phase = phase4
#     return result, phase, lc
