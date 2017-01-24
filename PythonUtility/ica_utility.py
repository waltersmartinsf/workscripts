
'''
Created September 28, 2016
Author: Walter Martins Filho

Goal: Help with the ICA, Independent Component Analysis
'''
import numpy as np
from sklearn.decomposition import PCA
import pandas as pd
import scipy #pearson correlation
import matplotlib.pyplot as plt

def rotate_axis(x):
    mean_value = np.mean(x)
    x = x - mean_value
    x_new = -1 * x#invert the sign
    return x_new + mean_value

def error_flux(hoststar_flux,hoststar_eflux,ref_star_flux,ref_star_eflux):
        _flux = hoststar_flux/ref_star_flux
        part1 = (hoststar_eflux/hoststar_flux)**2
        part2 = (ref_star_eflux/ref_star_flux)**2
        _eflux = _flux*np.sqrt(np.array(part1.values + part2.values))
        return error_flux

def pearson_ica_test(ica_signal,original_signal,save_dir):
    '''
    Return Pearson Statistcs about which column in the ica output is
    correlate with the main first principal component that corresponds
    to the light curve transit.

    Input:

    ica_signal: pandas dataframe
    original_signal: pandas dataframe

    '''

    pca = PCA(n_components=len(original_signal.columns))
    H = pca.fit_transform(original_signal)
    H = pd.DataFrame(H)

    H.plot(grid=True)
    print('Scatter 1st component = ',np.std(H[0]))
    plt.title('PCA Components')
    plt.savefig(save_dir+'PCA_components_.png')
    plt.savefig(save_dir+'PCA_components_.pdf')

    pearson,pvalue = np.zeros(ica_signal.shape[1]), np.zeros(ica_signal.shape[1])

    component_id = 0
    for i in range(ica_signal.shape[1]):
        pearson[i], pvalue[i] = scipy.stats.pearsonr(H[0],ica_signal[i])
        print(pearson[i], pvalue[i])
        if abs(pearson[i]) == abs(pearson).max():
            print('** Light curve on column = ',i,'\n')
            component_id = i
        else:
            print('** Probabily, this is not the light curve \n')
    return component_id
    #Test if we need to ratoate Y-axis on signal dataframe
#     print 'Testing the sign in the component ',component_id,' : \n'
#     pearson_ica, pvalue_ica = scipy.stats.pearsonr(original_signal[component_id],ica_signal[component_id].values)
#     if pearson_ica < 0:
#         print 'Rotate Y-axis:'
#         ica_signal[component_id] = rotate_axis(ica_signal[component_id])
