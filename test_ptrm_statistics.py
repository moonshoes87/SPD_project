#!/usr/bin/env python

import numpy
import unittest
import lib_ptrm_statistics as lib_ptrm



class CheckpTRMparams(unittest.TestCase):

    tmin = 20.
    tmax = 50.
    ptrm_x = [1., 2., 3., 5., 6.]
    ptrm_y = [7., 6., 4., 2.4, 0.]
    ptrm_temps = [10, 20, 30, 40, 50]
    ptrm_starting_temps = [20, 30, 40, 50, 60]
    ref_n = 3
    
    def test_get_n_ptrm(self):
        result = lib_ptrm.get_n_ptrm(self.tmin, self.tmax, self.ptrm_temps, self.ptrm_starting_temps)
        self.assertEqual(self.ref_n, result)


if __name__ == "__main__":
    unittest.main()
