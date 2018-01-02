import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from pandas import DataFrame


def k2_dftest(sample,savedir,title='KS Test',xlabel='x',ylabel='y',show=True,saveplot=True,name_plot='K2_test.pdf'):
    '''
    Apply the Kolmogorov-Sminorv Test between columns of a pandas dataframe 
    '''

    k2_results = np.zeros((len(sample.columns),len(sample.columns)))
    print(k2_results)

    for i in range(len(sample.columns)):
        for j in range(len(sample.columns)):
            x = stats.ks_2samp(sample[sample.columns[i]].values,sample[sample.columns[j]].values)
            print(i,j,x)
            k2_results[i][j] = 1- x.pvalue

    plt.figure(figsize=(14,9))
    plt.imshow(k2_results,cmap=plt.cm.coolwarm)
    for i in range(k2_results.shape[0]):
        for j in range(k2_results.shape[1]):
            plt.text(i,j,str(round(k2_results[i][j],6)),horizontalalignment='center',verticalalignment='center',color='white')
    
    plt.colorbar(label='p-value')
    names = sample.columns

    plt.xticks(np.arange(0,len(names),1),names)
    plt.yticks(np.arange(0,len(names),1),names,rotation=90,verticalalignment='center')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if saveplot == True:
        plt.savefig(savedir+name_plot)
    if show == True:
        plt.show()
    else:
        plt.close()

    return k2_results