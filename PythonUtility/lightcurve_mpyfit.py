'''
Author: Walter S. Martins-Filho
email:	walter @on.br @lpl.arizona.edu
		waltersmartinsf @gmail.com
		walterwsmf @outlook.com

This code fitting a lightcurve model to an exoplanetary transit.
It was based in Mandel and Agol (2002) for the mathematical
theory behind the light curve model.
'''
import os
import sys
import glob
import time

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm

import numpy as np
from scipy.stats import chisquare
from pandas import DataFrame, read_csv, concat

import emcee
import mpyfit
import occultation_bic

def lightcurve_fit(_time,_flux,_eflux,_am,_RpRs,_aR,_i,_u1,_u2,_P,_e,_omega):
	'''
	'''
	########################################################################################
	########################################################################################
	#Importing data to make the model
	Rp_Rs 			= float(_RpRs)
	limbB1, limbB2 	= float(_u1),float(_u2)
	inc 			= np.radians(float(_i))
	a_Rs 			= float(_aR)
	hjd 			= _time
	rawflux 		= _flux
	eflux 			= _eflux
	am 				= _am
	period 			= float(_P)
	omega 			= float(_omega)
	e 				= float(_e)

	########################################################################################
	########################################################################################
	#Functions based on Mandel and Agol (2002)
	def simple_model(JD, startpar0, startpar1):
		phase = (JD - startpar1)/period
		distance_vector = occultation_bic.delta(phase,inc) * a_Rs
		#model = occultation_fn(distance_vector,startpar0,limbB1,limbB2,show=False)
		model = occultation_bic.occultquad(distance_vector, limbB1, limbB2, startpar0, len(hjd))
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
    ########################################################################################
    ########################################################################################
    #Creating dictionary with guess parameters
	startpar = [float(Rp_Rs),np.mean(hjd), 1., 0.]
	print 'Parameters for fitting a phase model = ',startpar,'\n'

	PARINFO = [{'value':Rp_Rs,'limits':(0,1.)},{'value':np.mean(hjd)},{'value':1.},{'value':0.,'fixed':True}]
	pfit1, results1 = mpyfit.fit(residuals_am_exp, startpar, args = (hjd,rawflux,eflux), parinfo=PARINFO)
	#
	# model1 = model_am_exp(hjd,pfit1[0],pfit1[1],pfit1[2],pfit1[3])
	# phase1 = (hjd - pfit1[1])/period
	# #results = [{'bestfit':pfit1},{'fit_report':results1},{'phase':phase1},{'model':model1}]
	return pfit1
