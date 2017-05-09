import numpy as np
import ctypes
import requests


# define C double array
array_1d_double = np.ctypeslib.ndpointer(dtype=ctypes.c_double,ndim=1,flags=['C_CONTIGUOUS','aligned'])

# load library
lib_trans = np.ctypeslib.load_library('lib_transit.so','/Users/walter/github/workscripts/PythonUtility/lightcurveModel')

# load function occultquadC from transit library
occultquadC = lib_trans.occultquad
occultquadC.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, array_1d_double]
occultquadC.restype = None


if __name__ == "__main__":



    Rp = 0.055
    tmid = 0.9721 
    aR = 14.07
    i = 88.75
    u1 = 0.016
    u2 = 0.5029
    P = 3.336817
    e = 0
    omega = 0
    t = np.linspace(0.85,1.05,200)
    n = len(t)
    t = np.require(t,dtype=ctypes.c_double,requirements='C')
    curve = np.zeros(n,dtype=ctypes.c_double)
    occultquadC(t,Rp,aR,P,i,u1,u2,e,omega,tmid,n,curve) # old results


    import matplotlib.pyplot as plt
    plt.plot(t,curve)
    plt.show()
