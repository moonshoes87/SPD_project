#! /usr/bin/env python

import numpy
from numpy import linalg
#import copy
import spd



def get_zdata_segment(zdata, t_Arai, tmin, tmax):  # oh, right, change this to use start and end.  
    lower_bound = t_Arai.index(tmin)
    upper_bound = t_Arai.index(tmax)
    zdata_segment = zdata[lower_bound:upper_bound + 1] # inclusive of upper bound
    return zdata_segment

#ex_zdata = [[1., 1., 1.],[1., 2., 3.], [2., 4., 4.5], [2., 3., 7.5]]

def get_cart_primes_and_means(zdata_segment):
# seek simpler syntax
    zdata = numpy.array(zdata_segment)
    X1 = zdata[:,0]
    X2 = zdata[:,1]
    X3 = zdata[:,2]
    means = list(numpy.mean(zdata.T,axis=1))
    X1_prime = list(X1 - means[0])#X1_mean)
    X2_prime = list(X2 - means[1])
    X3_prime = list(X3 - means[2])
    return {'X1_prime': X1_prime, 'X2_prime': X2_prime, 'X3_prime': X3_prime, 'X1': X1, 'X2': X2, 'X3': X3, 'X1_mean': means[0], 'X2_mean': means[1], 'X3_mean': means[2], 'means': means }

def get_orientation_tensor(X1_p, X2_p, X3_p):
    X1_p, X2_p, X3_p = numpy.array(X1_p), numpy.array(X2_p), numpy.array(X3_p)
    orient_tensor = [[sum(X1_p * X1_p), sum(X1_p * X2_p), sum(X1_p * X3_p)],
                     [sum(X1_p * X2_p), sum(X2_p * X2_p), sum(X2_p * X3_p)],
                     [sum(X1_p * X3_p), sum(X2_p * X3_p), sum(X3_p * X3_p)]]
    tau, V = numpy.linalg.eig(orient_tensor)
    return {'orient_tensor': orient_tensor, 'tau': tau, 'V': V}

def get_eigenvectors_and_eigenvalues(orient_tensor):
    # may need to sort eigenvectors, eigenvalues.  possibly not in order...?
    pass

def get_reference_vector(X1_prime, X2_prime, X3_prime):
    n = len(X1_prime) - 1
    X1 = X1_prime[0] - X1_prime[n]
    X2 = X2_prime[0] - X2_prime[n]
    X3 = X3_prime[0] - X3_prime[n]
    print numpy.array([X1, X2, X3])
    return numpy.array([X1, X2, X3])
#    return numpy.array([0, 0, 0])
    

#M = (zdata_segment-mean(zdata_segment.T,axis=1)).T

