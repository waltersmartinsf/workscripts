def trunc(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    f = round(f, n)
    slen = len('%.*f' % (n, f))
    return str(f)[:slen] ###perottoni

def import_file_(file,columns):
    '''
    This function will import your NPY file and will return in a dataframe.

    INPUT:

    file: the name of the file, or the path+name of the file
    columns: the name of the columns in the file in a numpy.ndarray or a list

    OUTPUT:

    dataframe with the information in the file with the columns's name exactly how you set up.
    '''
    from pandas import DataFrame
    from numpy import load

    x = load(file)
    x = DataFrame(x)
    x.columns = columns
    return x

metals = [0.01,0.02]
age_values = [9.0]
name_file = 'xMIST_iso_metal_0.01_age_9.0.txt.npy'
name_columns = ['log10_isochrone','age','initial_mass']

x = import_file_(file=name_file,columns=name_columns)

print(x)

#plot example:
import matplotlib.pyplot as plt

plt.figure()
x.initial_mass.plot()

plt.figure()
x.log10_isochrone.plot()

plt.show()