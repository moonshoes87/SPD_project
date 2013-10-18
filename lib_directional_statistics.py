#! /usr/bin/env python

import numpy
#import copy
import spd
#import known_values
#import lib_arai_plot_statistics as lib_arai



def get_partial_zdata(zdata, t_Arai, tmin, tmax):
    lower_bound = t_Arai.index(tmin)
    upper_bound = t_Arai.index(tmax)
    partial_zdata = zdata[lower_bound:upper_bound + 1] # inclusive of upper bound
    return partial_zdata

def get_cart_averages(zdata):
    X = numpy.sum(zdata, axis=0)
    print X
    print len(zdata)
    X /= float(len(zdata))
    print X
    X1 = X[0]
    X2 = X[1]
    X3 = X[2]
    return {'X1': X1, 'X2': X2, 'X3': X3}



zdata = [[1., 1., 1.],[1., 2., 3.], [2., 4., 4.5], [2., 3., 7.5]]

def get_cart_primes(zdata = zdata, X1_avg = 1.5, X2_avg = 2.5, X3_avg = 4.):
    zdata = numpy.array(zdata)
    X1 = zdata[:,0]
    X2 = zdata[:,1]
    X3 = zdata[:,2]
    X1_prime = list(X1 - X1_avg)
    X2_prime = list(X2 - X2_avg)
    X3_prime = list(X3 - X3_avg)
    print zdata
    return {'X1_prime': X1_prime, 'X2_prime': X2_prime, 'X3_prime': X3_prime, 'X1': X1, 'X2': X2, 'X3': X3}

def get_orientation_tensor():
    pass
