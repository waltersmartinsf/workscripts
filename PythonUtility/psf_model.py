import numpy as np 
import matplotlib.pyplot as plt

'''
This code create a synthetic PSF in a 2D array.
'''

def makeGaussian(amplitude,x0=256,y0=512,sigx=4,sigy=4,CCDsize=[1024,1024],b=0):
    '''
    '''
    x = np.arange(0,CCDsize[0],1)
    y = np.arange(0,CCDsize[1],1)
    x, y = np.meshgrid(x,y)
    GX = (x-x0)**2 / (2*sigx**2)
    GY = (y-y0)**2 / (2*sigy**2)
    Gaussian2D = amplitude * np.exp(-(GX+GY)) + b
    return Gaussian2D
    
#Example.
g = makeGaussian(100)
plt.figure()
plt.imshow(g,origin='lower')
plt.colorbar()
plt.ylim(512-10,512+10)
plt.xlim(256-10,256+10)
plt.show()