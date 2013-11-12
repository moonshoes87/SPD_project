#!/usr/bin/env python

import numpy
import unittest
import lib_ptrm_statistics as lib_ptrm



class CheckpTRMparams(unittest.TestCase):

    tmin = 20.
    tmax = 50.
    x_ptrm = [1., 2., 3., 5., 6.]
    y_ptrm = [7., 6., 4., 2.4, 0.]
    ptrm_temps = [10, 20, 30, 40, 50]
    ptrm_starting_temps = [20, 30, 40, 50, 60]
    ref_n = 4.
    ref_steps = [10, 20, 30, 40]

    x_Arai = [1., 2., 2.5, 3., 4., 5., 6., 8.] # ptrm initially acquired at a temp
    t_Arai = [0,  10, 20,  30, 40, 50, 60, 70]
    ref_max_ptrm_check = 1.25
    ref_max_percent = 125.

    x_int = 8.5
    ref_delta_CK = 1.25 / 8.5 * 100
    
    def test_n_ptrm(self):
        result = lib_ptrm.get_n_ptrm(self.tmin, self.tmax, self.ptrm_temps, self.ptrm_starting_temps)
        self.assertEqual(self.ref_n, result[0])
        self.assertEqual(self.ref_steps, result[1])

    def test_max_ptrm_check(self):
#        [20, 30, 40] --> [2, 3, 5] # ptrm
#        [20, 30, 40] --> [2.5, 3., 4.]         # regular
        result = lib_ptrm.get_max_ptrm_check(self.ref_steps, self.ptrm_temps, self.x_ptrm, self.t_Arai, self.x_Arai)
        max_ptrm_diff, max_diff_percent = result[0], result[1]
        self.assertAlmostEqual(self.ref_max_ptrm_check, max_ptrm_diff)
        self.assertAlmostEqual(self.ref_max_percent, max_diff_percent)


    def test_delta_CK(self):
        result = lib_ptrm.get_delta_CK(self.ref_max_ptrm_check, self.x_int)
        self.assertAlmostEqual(self.ref_delta_CK, result)
#        ref_max_check  / x_int  * 100
        
    def test_DRAT(self):
        delta_y_prime = 5.
        delta_x_prime = 8.
        L = numpy.sqrt(25 + 64)
        ref_DRAT = self.ref_max_ptrm_check / L * 100.
        result = lib_ptrm.get_DRAT(delta_y_prime, delta_x_prime, self.ref_max_ptrm_check)
        self.assertAlmostEqual(ref_DRAT, result)


if __name__ == "__main__":
    unittest.main()
