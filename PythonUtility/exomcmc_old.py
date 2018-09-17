import numpy as np
from pandas import DataFrame
import emcee
import lightcurve #Kyle's lightcurve code
import lightcurveMCMC

def exomcmc(hjd,rawflux,eflux,airmass,pfit,sigma,ndim,nwalkers,Nsize,a_Rs,inc,period,ecc,omega,limbB1,limbB2,show=True,space=1e-4,lnf0=np.log(0.43)):
    '''
    Executing the Monte Carlo Markov Chain in Exoplanetary Transits
    ---
    INPUT:
    pfit:
    ndim
    nwalkers
    lnprob: A function that takes a vector in the parameter space as input and
            returns the natural logarithm of the posterior probability for that
            position
    show:   True or False. Shows in terminal the output in each step.
    space:  the size of the possible random value pick in the distribution
            on the space parameter. This give to you how much the random walk
            will be far to relative in the initial value gived.
    Nsize:  number of total chains in the MCMC Algoritm.
    '''
    # Lightcurve model based on occultquad routine in C/Python from Kreidberg
    from occultquad import occultquad
    import mpyfit
    am = np.array(airmass)

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

    def simple_model(JD, startpar0, startpar1):
        phase = (JD - startpar1)/period
        distance_vector = delta(phase,inc) * a_Rs
        # print distance_vector
        #model = occultation_fn(distance_vector,startpar0,limbB1,limbB2,show=False)
        model = occultquad(distance_vector, limbB1, limbB2, startpar0)#, len(hjd))
        model = model[0]
        return model

    def model_am_exp(hjd,startpar0,startpar1,startpar2,startpar3):
        model = simple_model(hjd,startpar0, startpar1)
        model_am = model * startpar2 * np.exp(-1. * startpar3 * am) #multiply light curve by factor x exponential
        return model_am

    def residuals_am_exp(params,args): #residuals function for mpfit
        RpRs = params[0]
        Tc = params[1]
        mu1 = params[2]
        mu2 = params[3]
        hjd, data, eps_data = args
        model = model_am_exp(hjd,RpRs,Tc,mu1,mu2)
        return (data-model)/eps_data

    def model_am_linear(hjd,startpar0,startpar1,startpar2,startpar3):
        model = simple_model(hjd,startpar0, startpar1)
        model_am = model * (startpar2 * am + startpar3) #multiply light curve by factor x exponential
        return model_am

    def residuals_linear(params,args): #residuals function for mpfit
        RpRs = params[0]
        Tc = params[1]
        mu1 = params[2]
        mu2 = params[3]
        hjd, data, eps_data = args
        model = model_am_linear(hjd,RpRs,Tc,mu1,mu2)
        return (data-model)/eps_data

    def model_am_2deg(hjd,startpar0,startpar1,startpar2,startpar3,startpar4):
        model = simple_model(hjd,startpar0, startpar1)
        model_am = model * (startpar2 + startpar3*am + startpar4*am) #multiply light curve by factor x exponential
        return model_am

    def residuals_2deg_mpfit(params,args): #residuals function for mpfit
        RpRs = params[0]
        Tc = params[1]
        mu1 = params[2]
        mu2 = params[3]
        mu3 = params[4]
        hjd, data, eps_data = args
        model = model_am_2deg(hjd,RpRs,Tc,mu1,mu2,mu3)
        return (data-model)/eps_data

    def lnlike(theta,x,y,yerr):
        """
        THis function returns the maximum likelihood following the ideas of Hoog et al.(2010)
        to fitting a model to data using Markov Chain Monte Carlo Method using Metropole-Hastings Algoritm and Gibbs Sample.
        ___
        INPUT

        theta: turple with the parameters to fit and the confindence interval of fractional size f
        x, y: data to fit
        yerr: error of our data

        OUTPUT:

        natural logaritm of the likelihood function
        """
        p1, p2, p3, lnf = theta

        model = model_am_exp(x,p1, p2, p3, 0)
        # fit_result_lnlike = lightcurveMCMC.lightcurve_fit(float(p1),float(p2),float(p3),a_Rs,inc,limbB1,limbB2,period,e,omega,np.array(x),np.array(y),np.array(yerr))
        # model_standart = (fit_result_lnlike.final_curve - np.mean(fit_result_lnlike.final_curve))/np.std(fit_result_lnlike.final_curve)
        # model = model_standart * np.std(y) + np.mean(y)

        inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
        return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

    def lnprior(theta):
        """
        This function is where our sicence prior is established.
        ___
        INPUT

        theta: turple with the parameters to fit and the confindence interval of fractional size f

        OUTPUT

        The cts of the with prior information
        """
        p1, p2, p3, lnf = theta

        minRp, maxRp = pfit[0] - 7.0*sigma[0], pfit[0] + 7.0*sigma[0]
        minTc, maxTc = pfit[1] - 7.5*sigma[1], pfit[1] + 7.5*sigma[1]
        minA, maxA = pfit[2] - 7.5*sigma[2], pfit[2] + 7.5*sigma[2]

        if minRp < p1 < maxRp and  minTc < p2 < maxTc and minA < p3 < maxA and -10.0 < lnf < 1.0:
            return 0.0
        return -np.inf

    def lnprob(theta, x, y, yerr):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(theta, x, y, yerr)

    ############################################################################
    ############################################################################
    hjd,rawflux,eflux= np.array(hjd),np.array(rawflux),np.array(eflux)
    #defining the veriables
    p1, p2, p3 = pfit
    s1, s2, s3 = sigma
    lnf = lnf0
    minTc, maxTc = p2 - 2.5* s2, p2 + 2.5*s2
    if show == True:
        print(minTc, maxTc)
    startpar = [p1, p2, p3,lnf]

    if show == True:
        print('Initial parameters = ',startpar)
        print('Initial Dispersion of the initial parameters = ',sigma)
    #ndim, nwalkers = 4, 8
    pos = [startpar + space*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(hjd, rawflux, eflux))
    sampler.run_mcmc(pos, Nsize)
    if show == True:
        print('Acceptance fraction = ',np.mean(sampler.acceptance_fraction))
    #Create the dataframe OUTPUT
    #samples = sampler.chain[:, 50:, :].reshape((-1, ndim)) #test with only the first 50 chains
    data = DataFrame(sampler.chain.reshape((-1, ndim)))
    data.columns = ['RpRs','Tc','A','lnf']
    if show == True:
        print('The percentiles for Rp/Rs are: ')
        print('31% Rp/Rs = ',np.percentile(data.RpRs,31),' and the inferior error are = ', abs(np.percentile(data.RpRs,50)-np.percentile(data.RpRs,31)))
        print('50% Rp/Rs = ',np.percentile(data.RpRs,50))
        print('81% Rp/Rs = ',np.percentile(data.RpRs,81),' and the superior error are = ', abs(np.percentile(data.RpRs,50)-np.percentile(data.RpRs,81)))

        print('The percentiles for Tc are: ')
        print('31% Tc = ',np.percentile(data.Tc,31),' and the inferior error are = ', abs(np.percentile(data.Tc,50)-np.percentile(data.Tc,31)))
        print ('50% Tc = ',np.percentile(data.Tc,50))
        print('81% Tc = ',np.percentile(data.Tc,81),' and the superior error are = ', abs(np.percentile(data.Tc,50)-np.percentile(data.Tc,81)))

        print('The percentiles for A are: ')
        print('31% A = ',np.percentile(data.A,31),' and the inferior error are = ', abs(np.percentile(data.A,50)-np.percentile(data.A,31)))
        print('50% A = ',np.percentile(data.A,50))
        print('81% A = ',np.percentile(data.A,81),' and the superior error are = ', abs(np.percentile(data.A,50)-np.percentile(data.A,81)))

        print('The percentiles for lnf are: ')
        print('31% lnf = ',np.percentile(data.lnf,31),' and the inferior error are = ', abs(np.percentile(data.lnf,50)-np.percentile(data.lnf,31)))
        print('50% lnf = ',np.percentile(data.lnf,50))
        print('81% lnf = ',np.percentile(data.lnf,81),' and the superior error are = ', abs(np.percentile(data.lnf,50)-np.percentile(data.lnf,81)))

    return data

def mcmc(hjd,rawflux,eflux,pfit,sigma,N,a_Rs,inc,limbB1,limbB2,period,e,omega,show=True):
    '''
    Metropolis-Hastings MCMC
    '''
    from numba import jit
    def chisquare_dof(x,y,eps,pfit,a_Rs,inc,limbB1,limbB2,period,e,omega):
        """
        Return the chisquare data.

        """
        a,b,c = pfit
        #residuos = (y - model_am_exp(hjd,a,b,c,0))/eps
        fit_result_lnlike = lightcurveMCMC.lightcurve_fit(float(a),float(b),float(c),a_Rs,inc,limbB1,limbB2,period,e,omega,np.array(x),np.array(y),np.array(eps))
        model_standart = (fit_result_lnlike.final_curve - np.mean(fit_result_lnlike.final_curve))/np.std(fit_result_lnlike.final_curve)
        model = model_standart * np.std(y) + np.mean(y)
        residuos = (y- model)/eps
        chi2 = sum(residuos**2)
        return chi2,residuos

    def pick_mcmc(startpar):
        pick = np.zeros(len(startpar))
        indx = range(len(startpar))
        pick[int(np.random.choice(indx))] = 1.
        return pick

    @jit
    def mcmc(x,y,eps,sigma,pfit,togsig,N,a_Rs,inc,limbB1,limbB2,period,e,omega):
        param = np.zeros(len(pfit))
        param = pfit + np.random.normal(size=len(pfit)) * sigma * togsig
        #if param[0] < 0: #thi is a force in our mcmc routine to maintenece the first parameter inside the phase space
        #    while param[0] < 0:
        #        param = pfit + np.random.normal(size=len(pfit)) * sigma * togsig
        chi2, residuos = chisquare_dof(x,y,eps,param,a_Rs,inc,limbB1,limbB2,period,e,omega)

        result = []
        chi2_result = []
        accept = 0
        for i in range(int(N)):
            pick = pick_mcmc(pfit)
            dparam = np.random.choice(pick)* np.random.normal(size=len(pfit)) * sigma *togsig
            param_new = param + dparam
            chi2new, residuos = chisquare_dof(x,y,eps,param_new,a_Rs,inc,limbB1,limbB2,period,e,omega)
            prob = np.exp(0.5*(chi2-chi2new))
            a = min([prob,1])
            u = np.random.random()
            if u <=a:
                accept = accept + 1
                result.append(param_new)
                chi2_result.append(chi2new)
                param = param_new
                chi2 = chi2new
        return result,chi2_result,accept/N
    ##
    chi2, residuos = chisquare_dof(hjd,rawflux,eflux,pfit,a_Rs,inc,limbB1,limbB2,period,e,omega)

    if chi2/(len(rawflux)-3.) > 1:
        print('Your chi-squared value is more than 1. Inflating error bars.')
        eflux = eflux * np.sqrt(chi2/(len(rawflux)-3.))

    if show == True:
        print('Input data: ')
        print('Iter = ',N)
        print('Start parameters: ',np.array(pfit))
        print('errors: ', np.array(sigma))
    Result,chi2result,rate = mcmc(hjd,rawflux,eflux,sigma,pfit,[3.,2.5,3.5],N,a_Rs,inc,limbB1,limbB2,period,e,omega)
    return Result,chi2result,rate