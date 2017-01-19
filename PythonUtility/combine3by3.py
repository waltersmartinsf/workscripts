'''
Script to call IRAF and combine images 3 by 3

---
Requisites:
- PyRAF
- Numpy
- Astroconda Environment with IRAF legacy and Py.2.7
- login.cl file from IRAF at the same directory
---


'''
################################################################################
################################################################################

from pyraf import iraf
import glob
import string

#edit this term every time
label_science = 'wasp33b'

################################################################################
################################################################################
def split_list(alist, wanted_parts=1):
    '''
    Function that split a array in the number of groups, call wanted_parts,
    that we want.

    ---
    INPUT:

    alist: list of elements that we want to split_list
    wanted_parts: number of smaller groups

    ---

    Example:

    A = [0,1,2,3,4,5,6,7,8]

    print len(A)/3
    print split_list(A, wanted_parts=1)
    print split_list(A, wanted_parts=2)
    print split_list(A, wanted_parts=len(A)/3)

    Output of the Example:

    3
    [[0, 1, 2, 3, 4, 5, 6, 7, 8]]
    [[0, 1, 2, 3], [4, 5, 6, 7, 8]]
    [[0, 1, 2], [3, 4, 5], [6, 7, 8]]

    ---
    Source:
    http://stackoverflow.com/questions/752308/split-list-into-smaller-lists
    '''
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
             for i in range(wanted_parts) ]

################################################################################
################################################################################

science_images = glob.glob('AB'+label_science+'*.fits')

if len(science_images) > 500:
    if len(science_images) % 2 == 0:
        print '\nCombining the science images for 2 by 2: \n'
        print 'Number of final images = ',len(science_images)/2,'\n'
        print 'Please, use the script to combine 2 by 2.'
    if len(science_images) % 3 == 0:
        print '\nCombining the science images for 3 by 3: \n'
        print 'Number of final images = ',len(science_images)/3,'\n'
        list_for_combine = split_list(list(science_images),wanted_parts=len(science_images)/3)
        for i in range(len(list_for_combine)):
            ablist = string.join(list_for_combine[i],',')
            iraf.imcombine(ablist,'CAB'+label_science+str(i).zfill(4)+'.fits')
    else:
        print '\nCan not combine using 3 by 3.\n'
