import numpy as np
import ctypes
import requests

'''
Author: Walter Martins-Filho
email: walter at on.br
       walter at lpl.arizona.edu

This code create a lightcurve to fit with a exoplaentary transit.
Here I use Kyle's C routine for that, describe in

MandelTransits.c
compile

files, and create the library lib_transit.so.

This code only makes the import parameters more easily. 

'''

def lighthcurve(exoinfo):
	'''
	This function create a synthetic lightcurve to fit the primary or
	secondary exoplanetary transit, given information of the orbit, mass
	and radius of teh exoplanet and radius and mass of the start.
	It's necessary to define an array for time and give the midtransit time.

	INPUT:

	All information will give by a dictionary with the follow keywords.

	RpRs:   ratio between the radius of the exoplanet and its host star
	aR:     semi-major axis divided by radius of the host star
	i:      inclination in degrees
	u1:     linear limb-darkening coefficient
	u2:     quadratic limb-darkening coefficient
	period: period in days 
	ecc:    eccentricity
    omega:  longitude of pericentre
    time:   time-series in correct units (seconds, days,...)
	MidT:   midtransit time in correct units that you are using.

	OUTPUT:

	curve: synthetic lightcurve normalized in one.
	'''
	# define C double array
	array_1d_double = np.ctypeslib.ndpointer(dtype=ctypes.c_double,ndim=1,flags=['C_CONTIGUOUS','aligned'])
	# load library
	lib_trans = np.ctypeslib.load_library('lib_transit.so','/Users/walterwsmf/Documents/PythonUtility')
	# load function occultquadC from transit library
	occultquadC = lib_trans.occultquad
	occultquadC.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, array_1d_double]
	occultquadC.restype = None
    Rp = exoinfo.RpRs
   	tmid = exoinfo.MidT 
   	aR = exoinfo.aR
   	i = exoinfo.i
   	u1 = exoinfo.u1
   	u2 = exoinfo.u2
   	P = exoinfo.period
    e = exoinfo.ecc
    omega = exoinfo.omega
    t = exoinfo.time
    n = len(t)
    t = np.require(t,dtype=ctypes.c_double,requirements='C')
    curve = np.zeros(n,dtype=ctypes.c_double)
    occultquadC(t,Rp,aR,P,i,u1,u2,e,omega,tmid,n,curve) # old results
	return curve