import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


def corner_beautiful_plot(data,bestfit,split,bins=50,labels=None,interpolation='nearest',cmap=plt.cm.gray,show=True):
    """
    Create a croner plot from a pandas dataframe input.
    ___
    INPUT:

    data:           pandas dataframe with N columns
    bestfit:        best fit values to mark in the corner plot as  horizontal and vertical lines
                    with N-values in a numpy array
    split:          boolean array with N-elements. This set if you want only the decimal part of
                    the values. Because this plot is thinking in exoplanetary transits, this was
                    created to remove the int part of the Julian Date
    bins:           the number of bins that you whant to show at the 2D-histogram of the data at
                    the lower subplots of the corner plot. Default is 50 bins. 
    labels:         labels of each column of the corner plot. N-list of strings. Default is None,
                    and the code will use the names of the columns in the pandas dataframe input.
    interpolation:  string, interpolation of the 2D-histogram. Default is 'nearest'. 
                    Possible options are 'none', 'bilinear', 'bicubic', 'spline16', 'spline36', 
                    'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 
                    'bessel', 'mitchell', 'sinc','lanczos'.
    cmap:           `~matplotlib.colors.Colormap`, optional, default: 'hot'. If None, cmap to rc 
                    `image.cmap` value. `cmap` is ignored when `X` has RGB(A) information.
    show:           boolean value: True or False. Default is True. If default, then the function
                    will show information for each step.


    """

    #Setting the plot default parameters:
    def init_plotting2():
        plt.rcParams['figure.figsize'] = (14.0,14.0)
        plt.rcParams['font.size'] = 14
        #plt.rcParams['font.family'] = 'Times New Roman'
        plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
        plt.rcParams['axes.titlesize'] = 2*plt.rcParams['font.size']
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
    
    init_plotting2() #initiallizing the plot parameters
    
    #import the shape of the dataframe input to a variable to use
    #for obtain the data and set the size of the corner plot
    shape = data.shape #[0]rows, [1]columns 
    
    #If split have some value equal to True, than the exactly column
    #will have the int part removed.
    for k in range(shape[1]):
        if split[k] == True:
            data.iloc[:,k] = np.modf(data.iloc[:,k].values)[0]
    
    #creating the corner plot structure:
    #f: figure output
    #axarr: array-plot, matrix plot with the size of the columns of dataframe per
    #columns of the dataframe N X N plots
    f, axarr = plt.subplots(shape[1], shape[1])#, sharex=True)#, sharey=True)
    
    #if show == True, this routine will show the information for each step before starting
    #the porcedure to create the plot i X j.
    if show == True:
        print('Creating corner plot, shape = ',shape[1],' per ',shape[1])

    #Creating the plot i X j
    for i in range(shape[1]):
        for j in range(shape[1]):
            #information about what type of plot will be created
            #If i ==j, then will be created a histogram of the i-column at the 
            #pandas dataframe input.
            #If i !=j, then will be created a numpy.histogram2d variable, and with that
            #we will use the output as a image to be treated as plt.imshow function
            if show == True:
                print('Subplot = ',i,j)
                if i == j:
                    print('Histogram of ',i)
                else:
                    print('Density plot of ',i,' per ',j)
            # i == j, Diagonal plots: Hitograms plots from i-columns at dataframe input.
            if i == j:
                #remove histogram grid.
                axarr[i][j].grid() 
                #create the numpy.histogram variable
                #normed=True to match with the plt.imshow image
                #bins setting to be at the sqrt-scale, sqrt(len(i-column)), to show
                #appropriated scale of the data set. 
                H,bins_edges = np.histogram(data.iloc[:,i].values,normed=True,bins='sqrt')
                #plot the histogram at same conditions of the numpy.histogram H variable
                axarr[i][j].hist(data.iloc[:,i].values,normed=True,bins='sqrt')
                #plot a vertical line of the bestfit i-value to zero to maximum 
                #of the H-variable and set the x-limits and y-limits 
                axarr[i][j].vlines(bestfit[i],0,H.max(),color='red')
                axarr[i][j].set_ylim(0,H.max())
                axarr[i][j].set_xlim(data.iloc[:,i].min(),data.iloc[:,i].max())
                #setting the text with the value of i-bestfit rounded with 4-decimal numbers
                #color setting match with the vertical line
                axarr[i][j].text(bestfit[i],H.mean(),str(round(bestfit[i],4)),color='red')
                #remove y-axis to clean the corner plot with no-needed information
                axarr[i][j].get_yaxis().set_visible(False)
                #remove x-ticks if the plots are not max(j) per max(j).
                #Those plots represent histograms that are above others plots
                #Else, the x-ticks will be created, with the name of the column in the pandas
                #dataframe input or if labels are given, those will be used. 
                if labels == None:
                    if i == int(shape[1]-1):
                        axarr[i][j].set_xlabel(data.columns[i])
                else:
                    if i == int(shape[1]-1):
                        axarr[i][j].set_xlabel(labels[i])
                if i != int(shape[1]-1):
                    axarr[i][j].get_xaxis().set_visible(False)
                if i == int(shape[1]-1):
                    axarr[i][j].get_xaxis().set_visible(True)
                    axarr[i][j].locator_params(axis='x',nbins=6)
            #i != j plots: image plot from plt.imshow
            #We will map the numpy.histogram2d from the j-column (x-axis) per i-column (yaxis)
            #The change of j to x-axis and i to y-axis is because to mach the x-label of the
            #histograms in the diagonal plots with the images plots. So, this will set 
            #correctly the corner axis.
            else:
                #Create the numpy.histogram2d object (image) H from the data set j-column per
                #i-column, and set the x- and y-tiks labels range xedges and yedges, respectively.
                #We will use the number of bins given from the input parameters. The default is
                #bins = 50. 
                H, xedges, yedges = np.histogram2d(data.iloc[:,j].values,data.iloc[:,i].values,bins=bins)
                #Setting the minimum  and maximum  of each x- and y-ticks, rounded with 4-decimal.
                xmin, xmax = round(xedges.min(),4),round(xedges.max(),4)
                ymin, ymax = round(yedges.min(),4),round(yedges.max(),4)
                #remove the grid of the image.
                axarr[i][j].grid(b=False)
                #plotting the image H with cmap and interpolation given by input parameters.
                #the minimum and maximum of the map is set to be the difference between
                #the mean value of H with plus or the difference with standart deviation of H.
                #The aspect is set to be auto to adjust the box of the image with the x- and
                #y-ticks.
                axarr[i][j].imshow(H,cmap=cmap,origin='lower',
                                   extent=[xmin,xmax,ymin,ymax],aspect='auto',
                                   vmin=np.mean(H)-np.std(H),vmax=np.mean(H)+np.std(H),
                                   interpolation=interpolation)
                #Make the countour plot of H.
                axarr[i][j].contour(H,origin='lower',extent=[xmin,xmax,ymin,ymax])
                #Setting the box limits to match with the histograms plots at the diagonal 
                axarr[i][j].set_xlim(xmin,xmax)
                axarr[i][j].set_ylim(ymin,ymax)
                #plot the vertical and horizontal lines with the values of the bestfit i and j
                #values, with red color.
                axarr[i][j].hlines(bestfit[i],xmin,xmax,color='red')
                axarr[i][j].vlines(bestfit[j],ymin,ymax,color='red')
                #set the aspect of the box x- and y-ticks to auto too to match with the 
                #plt.imshow parameter.
                axarr[i][j].set_aspect('auto')
                #Fix the number of values at the array to clean the x- and y-ticks.
                #here, I choose to make only 6 bins, and this will give 5-ticks at maximum
                axarr[i][j].locator_params(axis='x',nbins=6)
                axarr[i][j].locator_params(axis='y',nbins=6)
                #Set the labels of each plot. If the labels are given as one of the input
                #parameters, these will be used instead of the names of the pandas 
                #dataframe input 
                if labels == None:
                    axarr[i][j].set_ylabel(data.columns[i])
                    axarr[i][j].set_xlabel(data.columns[j])
                else:
                    axarr[i][j].set_ylabel(labels[i])
                    axarr[i][j].set_xlabel(labels[j])
                #If the plot is between two other plots,
                #the y-ticks or the x-ticks will be removed.
                if (j > 0):
                    axarr[i][j].get_yaxis().set_visible(False)
                if (i < shape[1]-1):
                    axarr[i][j].get_xaxis().set_visible(False)
    #Remove the upper-plots of the diagonal plots
    #this create a triangle plot 
    for i, j in zip(*np.triu_indices_from(axarr, 1)): 
        axarr[i, j].set_visible(False)
    return f