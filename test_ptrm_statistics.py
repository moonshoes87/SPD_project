#!/usr/bin/env python

import numpy
import unittest
import spd
import lib_ptrm_statistics as lib_ptrm



class CheckpTRMparams(unittest.TestCase):

    tmin = 20.
    tmax = 50.
    x_ptrm = [1., 2., 3., 9., 6.]
    y_ptrm = [7., 6., 4., 2.4, 0.]
    ptrm_temps = [10, 20, 30, 40, 50]
    ptrm_starting_temps = [20, 30, 40, 50, 60]
    ref_n = 4.
    ref_steps = [10, 20, 30, 40]

    x_Arai = [1., 1., 2.5, 5.5, 7., 5., 6., 8.] # ptrm initially acquired at a temp
    t_Arai = [0,  10, 20,  30, 40, 50, 60, 70]
    #    1,   2,   3,  9.
    #  - 1,   2.5, 5.5, 7.
    #    0., -.5,  -2.5, 2.
    ref_max_ptrm_check = 2.5
    ref_sum_ptrm_check = abs(-1.)
    ref_check_percent = (2.5/ 5.5) * 100.
    ref_sum_abs_ptrm_check = 5.

    x_int = 8.5
    ref_delta_CK = 2.5 / 8.5 * 100

    delta_y_prime = 5.
    delta_x_prime = 8.

    ref_L = numpy.sqrt(25 + 64)


    
    def test_n_ptrm(self):
        result = lib_ptrm.get_n_ptrm(self.tmin, self.tmax, self.ptrm_temps, self.ptrm_starting_temps)
        self.assertEqual(self.ref_n, result[0])
        self.assertEqual(self.ref_steps, result[1])

    def test_max_ptrm_check(self):
#def get_max_ptrm_check(ptrm_checks_segment, ptrm_checks, ptrm_x, t_Arai, x_Arai):
        result = lib_ptrm.get_max_ptrm_check(self.ref_steps, self.ptrm_temps, self.x_ptrm, self.t_Arai, self.x_Arai)
        max_ptrm_diff, sum_ptrm_diffs, check_percent, sum_abs_ptrm_diffs = result[0], result[1], result[2], result[3]
        self.assertAlmostEqual(self.ref_max_ptrm_check, max_ptrm_diff)
        self.assertAlmostEqual(self.ref_sum_ptrm_check, sum_ptrm_diffs)
        self.assertAlmostEqual(self.ref_check_percent, check_percent)
        self.assertAlmostEqual(self.ref_sum_abs_ptrm_check, sum_abs_ptrm_diffs)


    def test_delta_CK(self):
        result = lib_ptrm.get_delta_CK(self.ref_max_ptrm_check, self.x_int)
        self.assertAlmostEqual(self.ref_delta_CK, result)

        
    def test_DRAT(self):
        ref_DRAT = (self.ref_max_ptrm_check / self.ref_L) * 100.
        DRAT, L = lib_ptrm.get_DRAT(self.delta_y_prime, self.delta_x_prime, self.ref_max_ptrm_check)
        self.assertAlmostEqual(ref_DRAT, DRAT)
        self.assertAlmostEqual(self.ref_L, L)

    def test_max_DEV(self):
        result = lib_ptrm.get_max_DEV(self.delta_x_prime, self.ref_max_ptrm_check)
        ref_max_DEV = (2.5 / 8.) * 100
        self.assertAlmostEqual(ref_max_DEV, result)
        
    def test_CDRAT(self):
        CDRAT, CDRAT_prime = lib_ptrm.get_CDRAT(self.ref_L, self.ref_sum_ptrm_check, self.ref_sum_abs_ptrm_check)
        ref_CDRAT, ref_CDRAT_prime = (1. / self.ref_L) * 100., (5. / self.ref_L) * 100
        self.assertAlmostEqual(ref_CDRAT, CDRAT)
        self.assertAlmostEqual(ref_CDRAT_prime, CDRAT_prime)
        
    def test_DRATS(self):
        #ref_DRATS = .9
        ref_DRATS = (self.ref_sum_ptrm_check / 7.) * 100.
        end = 4
        DRATS = lib_ptrm.get_DRATS(self.ref_sum_ptrm_check, self.x_Arai, end)
        self.assertAlmostEqual(ref_DRATS, DRATS)
        
    def test_DRATS_real_data(self):
        ref_drats = .9
        DRATS = spd.thing.get_DRATS()
        self.assertEqual(ref_drats, DRATS)
        


if __name__ == "__main__":
    unittest.main()
