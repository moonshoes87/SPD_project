#! /usr/bin/env python

from scipy import *
#import numpy



def get_vds(zdata, delta_y_prime, start, end):  # 
    """takes zdata array: [[1, 2, 3], [3, 4, 5]], delta_y_prime: 1, start value and end value.  gets vds and f_vds, etc. """
    print "calling get_vds()"
    vector_diffs = []
    for k in range(len(zdata)-1):
            # gets diff between two vectors                                                                      
        vector_diffs.append(sqrt(sum((array(zdata[k + 1 ])-array(zdata[k]))**2)))
    vector_diffs.append(sqrt(sum(array(zdata[-1])**2))) # last vector of the vds                             
    last_vector = sqrt(sum(array(zdata[-1])**2))
    vds = sum(vector_diffs)
    f_vds = abs(delta_y_prime / vds) # fvds varies, because of delta_y_prime, but vds does not.              
    vector_diffs_segment = vector_diffs[start:end]
    partial_vds = sum(vector_diffs_segment)
    max_diff = max(vector_diffs_segment)
    GAP_MAX = max_diff / partial_vds
    return {'max_diff': max_diff, 'vector_diffs': vector_diffs, 'specimen_vds': vds, 'f_vds': f_vds, 'vector_diffs_segment': vector_diffs_segment, 'partial_vds': partial_vds, 'GAP-MAX': GAP_MAX}




def get_FRAC(vds, vector_diffs_segment):   
    for num in vector_diffs_segment:
        if num < 0:
            raise ValueError
    FRAC=sum(vector_diffs_segment)/ vds
    print FRAC
    return FRAC

def get_R_corr2(x_avg, y_avg, x_segment, y_segment): # 
    numerator = 0
    denominator_x = 0
    denominator_y = 0
    for num, x in enumerate(x_segment):
        r = ((x - x_avg) **2 ) * ((y_segment[num] - y_avg) **2 )
        numerator += r
#    print "numerator", numerator
    for x in x_segment:
        denominator_x += ((x - x_avg) ** 2)
#    print "den_x", denominator_x
    for y in y_segment:
        denominator_y += ((y - y_avg) ** 2)
#    print "den_y", denominator_y
    R_corr2 = numerator / (denominator_x * denominator_y)
    print R_corr2
    return R_corr2

def get_R_det2(y_segment, y_avg, y_prime):
    """
    takes in an array of y values, the mean of those values, and the array of y prime values.  returns R_det2
    """
    top = 0
    for num, y in enumerate(y_segment):
        result = (y - y_prime[num]) ** 2
        top += result
    bottom = 0
    for num, y in enumerate(y_segment):
        result = (y - y_avg) **2
        bottom += result
    print "top, bottom", top, bottom
    R_det2 = 1 - (top / bottom)
    print R_det2
    return R_det2
