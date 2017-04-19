def lightcurve_fit(_RpRs,_tmid,A0,_aR,_i,_u1,_u2,_P,_e,_omega,_time,_flux,_eflux):
	'''
	Use lmfit package

	INPUT:
		_RpRs
	    _tmid: mid transit guess
	    _aR:
	    _i: inclination in degrees
	    _u1:
	    _u2:
	    _P:
	    _e:
	    _omega:
	    _time:
	    _flux:
	    _eflux:
	'''
	import numpy as np
	import ctypes
	import matplotlib.pyplot as plt
	from lmfit import minimize as lminimize
	from lmfit import Parameters, Parameter, report_fit
	from os import environ

	# LOAD IN C OCCULT FUNCTION
	array_1d_double = np.ctypeslib.ndpointer(dtype=ctypes.c_double,ndim=1,flags=['C_CONTIGUOUS','aligned'])

	# load library
	lib_test = np.ctypeslib.load_library("lib_transit.so",environ['UTILPATH'])

	# load fn from library and define inputs
	occultquadC = lib_test.occultquad
	occultquadC.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, array_1d_double]
	occultquadC.restype = None




	######################################################################
	#             BASIC TRANSIT FITTER (Rp,Tmid,A0,A1,A2) only
	######################################################################
	class param:
	    rp = 0
	    ar = 1
	    P = 2
	    i = 3
	    u1 = 4
	    u2 = 5
	    e = 6
	    w = 7
	    tmid = 8


	class lc_fitter(object):

	    def __init__(self,t,data,dataerr=None,init=None,bounds=None,errortype='chi2',airmass=None):

	        self.t = np.array(t)
	        self.y = np.array(data)

	        self.p_init = init
	        self.errortype = errortype
	        self.bounds = bounds
	        self.airmass = airmass

	        if type(dataerr) == type(None):
	            self.yerr = np.ones(len(t))
	        else:
	            self.yerr = dataerr

	        self.fit_lm()


	    def plot_results(self):

	        f,ax = plt.subplots(1)

	        ax.errorbar(self.t, self.y/self.amcurve,yerr=self.yerr,fmt='ko', alpha=.2, label="f + noise")
	        ax.plot(self.t, self.final_curve ,'--', c='r', lw=2., label="fitted f")

	        ax.set_xlabel("t")
	        ax.set_ylabel("f(t)")
	        ax.legend(loc="best")
	        #pl.xlim([-xlim*1.01,xlim*1.01])

	        plt.show()


	    def fit_lm(self):
	        # use Levenberg Mardquart method


	         # define objective function: returns the array to be minimized
	        def fcn2min(params, x, y, yerr):

	            n = len(x)
	            model = np.zeros(n,dtype=ctypes.c_double)
	            model = np.require(model,dtype=ctypes.c_double,requirements='C')

	            occultquadC( x,params['RpRs'].value,params['aRs'].value, params['period'].value, params['inc'].value,
	                        params['gamma1'].value, params['gamma2'].value, params['ecc'].value, params['omega'].value,
	                        params['tmid'].value, n, model )

	            model *= (params['a0'] + x*params['a1'] + x*x*params['a2'])

	            return (model - y)/yerr

	        #Rp,aR,P,i,u1,u2,e,omega,tmid,a0,a1,a2 = self.p_init
	        v = [ (i[0] != i[1]) for i in self.bounds ] # boolean array to vary parameters
	        pnames = ['RpRs','aRs','period','inc','gamma1','gamma2','ecc','omega','tmid','a0','a1','a2']
	        params = Parameters()

	        for j in range(len(self.p_init)):

	            # algorithm does not like distance between min and max to be zero
	            if v[j] == True:
	                params.add(pnames[j], value= self.p_init[j], vary=v[j], min=self.bounds[j][0], max=self.bounds[j][1] )
	            else:
	                if (self.bounds[j][0] == None):

	                    if (self.bounds[j][1] == None): # no upper bound
	                        params.add(pnames[j], value= self.p_init[j], vary=True )
	                    else: # upper bound
	                        params.add(pnames[j], value= self.p_init[j], vary=True,max = self.bounds[j][1] )

	                elif (self.bounds[j][1] == None):

	                    if (self.bounds[j][0] == None): # no lower bound
	                        params.add(pnames[j], value= self.p_init[j], vary=True )
	                    else: # lower bound
	                        params.add(pnames[j], value= self.p_init[j], vary=True,min = self.bounds[j][0] )

	                else:
	                    params.add(pnames[j], value= self.p_init[j], vary=v[j] )


	         # do fit, here with leastsq model
	        result = lminimize(fcn2min, params, args=(self.t,self.y,self.yerr))


	        params = result.params
	        n = len(self.t)
	        model = np.zeros(n,dtype=ctypes.c_double)
	        model = np.require(model,dtype=ctypes.c_double,requirements='C')
	        occultquadC( self.t,params['RpRs'].value,params['aRs'].value, params['period'].value, params['inc'].value,
	                    params['gamma1'].value, params['gamma2'].value, params['ecc'].value, params['omega'].value,
	                    params['tmid'].value, n, model )

	        self.final_model = model
	        self.residuals = result.residual
	        self.params = result.params
	        self.result = result


	        A0 = params['a0'].value
	        A1 = params['a1'].value
	        A2 = params['a2'].value
	        self.amcurve = A0 + self.t*A1 + self.t*self.t*A2
	        self.final_curve = self.final_model/self.amcurve
	        self.phase = (self.t-params['tmid'].value)/params['period']


	        # write error report
	        # report_fit(params)

	# if __name__ == "__main__":

	    # if running outside folder:
	    # from util import lc_fitter

	    ############################################################
	    ######################## Example  ##########################
	    # generate noisy data with known coefficients
	    # Rp = 0.055
	    # tmid = 0.9721 # mid transit guess
	    # aR = 14.07
	    # i = 88.75
	    # u1 = 0
	    # u2 = 0.5029
	    # P = 3.336817
	    # e = 0
	    # omega = 0
	    # a0 = 1
	    # a1 = np.random.normal(0.,0.001)
	    # a2 = np.random.normal(0.,0.001)
	    # t = np.linspace(0.7,1.2,200)


	    # n = len(t)
	    # t = np.require(t,dtype=ctypes.c_double,requirements='C')
	    # curve = np.zeros(n,dtype=ctypes.c_double)
	    # occultquadC(t,Rp,aR,P,i,u1,u2,e,omega,tmid,n,curve) # old results
	    # curve = curve * (a0 + a1*t + a2*t*t)
	    # data = curve + np.random.normal(0, 5e-4, len(curve))
	    # dataerr = abs( 5e-4+np.random.normal(0,5e-4,len(curve)) )
	##############################################################
	##############################################################
	######################## Input data  #########################
	t 		= _time
	Rp 		= _RpRs
	tmid 	= _tmid
	aR 		= _aR
	i 		= _i
	u1 		= _u1
	u2 		= _u2
	P 		= _P
	e 		= _e
	omega 	= _omega
	data 	= _flux
	dataerr = _eflux
	##############################################################
	######################## Fitting  ############################

	p_init = [Rp,aR,P,i,u1,u2,e,omega,tmid,A0,0,0]
	mybounds = [ (0.0,1), # Rp/Rs
	            (aR,aR), # aR
	            (P,P), # P
	            (i,i), # i
	            (u1,u1), # u1
	            (u2,u2), # u2
	            (e,e), # e
	            (omega,omega), # omega
	            (min(t),max(t)), # tmid
	            (None,None), # A0  #cooeficient linear for airmass correction
	            (0,0), # A1 #coefficient quadratic for aimass correction
	            (0,0)] # A2 #coefficient polynomial for airmass correction


	myfit = lc_fitter(t,data,dataerr=dataerr,init= p_init,bounds= mybounds,)
	# print(myfit.residuals)
	# myfit.plot_results()
	# report_fit(myfit.params)
	# print('my best Rp/Rs =',myfit.params['RpRs'],'+-',myfit.params['RpRs'].stderr)
	return myfit#.params,model
