#! /usr/bin/env python

import numpy
#import copy
import spd



def get_partial_zdata(zdata, t_Arai, tmin, tmax):  # oh, right, change this to use start and end.  
    lower_bound = t_Arai.index(tmin)
    upper_bound = t_Arai.index(tmax)
    partial_zdata = zdata[lower_bound:upper_bound + 1] # inclusive of upper bound
    return partial_zdata

#ex_zdata = [[1., 1., 1.],[1., 2., 3.], [2., 4., 4.5], [2., 3., 7.5]]

def get_cart_primes_and_means(zdata):
# seek simpler syntax
    zdata = numpy.array(zdata)
    X1 = zdata[:,0]
    X2 = zdata[:,1]
    X3 = zdata[:,2]
    means = list(numpy.mean(zdata.T,axis=1))
    X1_prime = list(X1 - means[0])#X1_mean)
    X2_prime = list(X2 - means[1])
    X3_prime = list(X3 - means[2])
    return {'X1_prime': X1_prime, 'X2_prime': X2_prime, 'X3_prime': X3_prime, 'X1': X1, 'X2': X2, 'X3': X3, 'X1_mean': means[0], 'X2_mean': means[1], 'X3_mean': means[2], 'means': means }

def get_orientation_tensor(zdata, means):
    
    pass

#M = (zdata_segment-mean(zdata_segment.T,axis=1)).T

