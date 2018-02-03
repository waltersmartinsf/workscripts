import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

def init_plotting(x=18,y=14):
    plt.rcParams['figure.figsize'] = (x,y)
    plt.rcParams['font.size'] = 20
    #plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 0.75*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = 0.65*plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.loc'] = 'best'
    plt.rcParams['axes.linewidth'] = 1

#init_plotting()


def plot_residuals(x,y1,y2,y1err=None,title='Plot',y1label='y1',xlabel='x',y2label='y2'):
    f = plt.figure()
    plt.suptitle(title)
    gs1 = GridSpec(2, 2, width_ratios=[1,2],height_ratios=[4,1])
    gs1.update(wspace=0.5)
    ax1 = plt.subplot(gs1[:-1, :])
    ax2 = plt.subplot(gs1[-1, :])
    ax1.grid()
    ax1.errorbar(x,y1,yerr=y1err,ecolor='g')
    ax1.set_xticklabels([])
    ax1.set_ylabel(y1label)
    ax1.set_xlim(x.min(),x.max())
    ax2.grid()
    ax2.plot(x,y2,color='green')
    plt.yticks(np.array([round(y2.min(),3), round((y2.min()+y2.max())/2.,3), round(y2.max(),3)]))
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(y2label)
    ax2.set_xlim(x.min(),x.max())
    plt.close()
    return f

def A4init_plotting():
    plt.rcParams['figure.figsize'] = (11.69,8.27)
    plt.rcParams['font.size'] = 15
    #plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 0.75*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = 0.65*plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.loc'] = 'best'
    plt.rcParams['axes.linewidth'] = 1