import numpy as np


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

def mcmc():
    