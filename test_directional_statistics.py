#! /usr/bin/env python

import unittest
import numpy
#import copy
import spd
#import known_values
import lib_directional_statistics as lib_direct


thing = spd.thing

class CheckDecInc(unittest.TestCase):
    
    zdata = [[0., 0., 0.],[1., 1., 1.],[1., 2., 3.], [2., 4., 4.5], [2., 3., 7.5], [0., 0., 0.]]
    partial_zdata = [[1., 1., 1.],[1., 2., 3.], [2., 4., 4.5], [2., 3., 7.5]]
    t_Arai = [1.,1.5, 3, 8, 10, 11]
    tmin = 1.5
    tmax = 10.

    X1_avg = 1.5
    X2_avg = 2.5
    X3_avg = 4.

    X1_prime = [-.5, -.5, .5, .5]
    X2_prime = [-1.5, -.5, 1.5, .5]
    X3_prime = [-3., -1., .5, 3.5]

    X1 = [ 1.,  1.,  2.,  2.] 
    X2 = [ 1.,  2.,  4.,  3.] 
    X3 = [ 1.,   3.,   4.5,  7.5]


    X1_p = numpy.array(X1_prime)
    X2_p = numpy.array(X2_prime)
    X3_p = numpy.array(X3_prime)
    # probably can write this more succintly, but good visualization
    orient_tensor = [[sum(X1_p * X1_p), sum(X1_p * X2_p), sum(X1_p * X3_p)],
                     [sum(X1_p * X2_p), sum(X2_p * X2_p), sum(X2_p * X3_p)],
                     [sum(X1_p * X3_p), sum(X2_p * X3_p), sum(X3_p * X3_p)]]
            

    def test_partial_zdata(self):
        result = lib_direct.get_partial_zdata(self.zdata, self.t_Arai, self.tmin, self.tmax)
        self.assertEqual(result, self.partial_zdata)

    def test_cart_averages(self):
        result = lib_direct.get_cart_primes_and_means(self.partial_zdata)
        self.assertEqual(result['X1_mean'], self.X1_avg)
        self.assertEqual(result['X2_mean'], self.X2_avg)
        self.assertEqual(result['X3_mean'], self.X3_avg)

    def test_cart_primes(self):
        result = lib_direct.get_cart_primes_and_means(self.partial_zdata)
        self.assertEqual(result['X1_prime'], self.X1_prime)
        self.assertEqual(result['X2_prime'], self.X2_prime)
        self.assertEqual(result['X3_prime'], self.X3_prime)

    def test_orientation_tensor(self):
        pass
       # result = lib_direct.get_orientation_tensor(self.X1_prime, self.X2_prime, self.X3_prime)
        
        

# see L5772 in thellier_gui_spd


#>>> for num, x in enumerate(X2):
#...     t = x * X2_prime[num]
#...     total += t



if __name__ == "__main__":
    unittest.main()
