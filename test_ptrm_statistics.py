#!/usr/bin/env python

import numpy
import unittest
import spd
import lib_ptrm_statistics as lib_ptrm
import lib_directional_statistics as lib_direct



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
        ref_DRATS = (1. / 7.) * 100.
        ref_DRATS_prime = (5. / 7.) * 100.
        end = 4
        DRATS, DRATS_prime = lib_ptrm.get_DRATS(self.ref_sum_ptrm_check, self.ref_sum_abs_ptrm_check, self.x_Arai, end)
        self.assertAlmostEqual(ref_DRATS, DRATS)
        self.assertAlmostEqual(ref_DRATS_prime, DRATS_prime)
        
    def test_DRATS_real_data(self):
        ref_drats = 0.928840447566
        DRATS = spd.thing.get_DRATS()
        self.assertAlmostEqual(ref_drats, DRATS)

    def test_mean_DRAT(self):
#        (1 / nPTRM) * (ref_sum_ptrm_check / L) = 1/4 * 1/ numpy.sqrt(25 + 64)
#        (1 / nPTRM) * (ref_abs_sum_ptrm_check / L)
#ref_L = numpy.sqrt(25 + 64)        
#    ref_sum_ptrm_check = abs(-1.)
#    ref_sum_abs_ptrm_check = 5.
#  ref_n = 4
        ref_mean_DRAT = 0.026499947000159001 * 100.
        ref_mean_DRAT_prime = 0.13249973500079501 * 100.
        mean_DRAT, mean_DRAT_prime = lib_ptrm.get_mean_DRAT(self.ref_sum_ptrm_check, self.ref_sum_abs_ptrm_check, self.ref_n, self.ref_L)
        self.assertAlmostEqual(ref_mean_DRAT, mean_DRAT)
        self.assertAlmostEqual(ref_mean_DRAT_prime, mean_DRAT_prime)

    def test_mean_DEV(self):
        ref_mean_DEV = (1. / 4.) * ( 1.  / 8.)  * 100
        ref_mean_DEV_prime = (1./ 4.) * (5. / 8.)  * 100
        mean_DEV, mean_DEV_prime = lib_ptrm.get_mean_DEV(self.ref_sum_ptrm_check, self.ref_sum_abs_ptrm_check, self.ref_n, self.delta_x_prime)
        self.assertAlmostEqual(ref_mean_DEV, mean_DEV)
        self.assertAlmostEqual(ref_mean_DEV_prime, mean_DEV_prime)

class CheckDeltaPal(unittest.TestCase):

    
    # add in initial TRM [0, 79., 11., 0, 1]
    PTRMS = numpy.array([[0, 79., 11., 0., 1.],[10, 78.05582281,  10.65530605,   0, 1], [20, 78., 10., 0., 0], [30, 60.5021261 , 8.24033996, 0, 1],[40, 65.37643521,  10.72414794,   0, 0.], [50., 70.60218755,   7.5674351, 0, 1. ]])

    TRM_1 = numpy.array([ 0.1873033 ,  0.96359193,  0.190809  ])

# cartesian  # and not including TRM_1
    PTRMS_cart = numpy.array([[ 0.20339007,  0.96148034,  0.18490007],[ 0.20475305,  0.96328734,  0.17364818], [ 0.4873076 ,  0.86138785,  0.14332577], [ 0.40937761,  0.89318751,  0.18608073], [ 0.32923249,  0.93502028,  0.131693]])

# direction
    PTRMS_dir = numpy.array([[ 78.05582281,  10.65530605,   1.], [78., 10., 1], [ 60.5021261 ,   8.24033996,   1.], [ 65.37643521,  10.72414794,   1.], [ 70.60218755,   7.5674351 ,   1. ]])    

# format of PTRMS = [[373.0, 144.66319008038636, -22.857955141406837, 6.5034930359418245e-11, 0], ...
# temp, dec, inc, moment, ZI or IZ

    PTRM_Checks_cart = numpy.array([[ 0.19245009,  0.96225045,  0.19245009], [ 0.44232587,  0.88465174,  0.14744196],[ 0.43643578,  0.87287156,  0.21821789],[ 0.27216553,  0.95257934,  0.13608276]])

    PTRM_Checks_dir = numpy.array([[ 78.69006753,  11.09580328,   1.], [ 63.43494882,   8.47871315,   1.], [ 63.43494882,  12.60438265,   1.], [ 74.0546041 ,   7.82123551,   1.]])

    PTRM_Checks = numpy.array([[10, 78.69006753,  11.09580328, 0],[30, 63.43494882,   8.47871315,   0], [40, 63.43494882,  12.60438265,  0], [50, 74.0546041 ,   7.82123551,   0]])

# format of PTRM_Checks
# temp, dec, inc, moment
#[[373.0, 276.84693928434058, 83.278957773601846, 6.0140429219203115e-11], ....


    correct_diffs = numpy.array([[  1.09399800e-02,  -7.70110000e-04,  -7.55002000e-03],[0,0,0], [  4.49817300e-02,  -2.32638900e-02,  -4.11619000e-03], [ -2.70581700e-02,   2.03159500e-02,  -3.21371600e-02], [  5.70669600e-02,  -1.75590600e-02,  -4.38976000e-03]])


    correct_C = numpy.array([[ 0.01093998, -0.00077011, -0.00755002],[ 0.01093998, -0.00077011, -0.00755002], [ 0.05592171, -0.024034  , -0.01166621], [ 0.02886354, -0.00371805, -0.04380337], [ 0.0859305 , -0.02127711, -0.04819313]])

    TRM_star = numpy.array([[ 0.1873033 ,  0.96359193,  0.190809  ],[ 0.21433005,  0.96071023,  0.17735005], [ 0.21569303,  0.96251723,  0.16609816], [ 0.54322931,  0.83735385,  0.13165956], [ 0.43824115,  0.88946946,  0.14227736],  [ 0.41516299,  0.91374317,  0.08349987]])


    def test_delta_pal_PTRM_vectors(self):
        ptrms_vectors, ptrms_checks_vectors, TRM_1 = lib_ptrm.get_delta_pal_vectors(self.PTRMS, self.PTRM_Checks)
        for num, vector in enumerate(ptrms_vectors):
            for n, i in enumerate(vector):
                print i, self.PTRMS_cart[num][n]
                self.assertAlmostEqual(i, self.PTRMS_cart[num][n])

    def test_delta_pal_TRM_1(self):
        ptrms_vectors, ptrms_checks_vectors, TRM_1 = lib_ptrm.get_delta_pal_vectors(self.PTRMS, self.PTRM_Checks)
        for num, value in enumerate(TRM_1.ravel()): # flattens array to avoid unnecessary looping
            self.assertAlmostEqual(self.TRM_1[num], value)
    
    def test_delta_pal_check_vectors(self):
        ptrms_vectors, ptrm_checks_vectors, TRM_1 = lib_ptrm.get_delta_pal_vectors(self.PTRMS, self.PTRM_Checks)
        for num, vector in enumerate(ptrm_checks_vectors):
            for n, i in enumerate(vector):
                print i, self.PTRM_Checks_cart[num][n]
                self.assertAlmostEqual(i, self.PTRM_Checks_cart[num][n])


    def test_diffs(self):
        diffs, C = lib_ptrm.get_diffs(self.PTRMS_cart, self.PTRM_Checks_cart, self.PTRMS, self.PTRM_Checks)
        for num, diff in enumerate(diffs):
            for n, v in enumerate(diff):
#                print self.correct_diffs[num][n], v
                self.assertAlmostEqual(self.correct_diffs[num][n], v)

    def test_C(self):
        diffs, C = lib_ptrm.get_diffs(self.PTRMS_cart, self.PTRM_Checks_cart, self.PTRMS, self.PTRM_Checks)
        for num, total in enumerate(C):
            for n, v in enumerate(total):
                self.assertAlmostEqual(self.correct_C[num][n], v)

    def test_TRM_star(self):
        TRM_star = lib_ptrm.get_TRM_star(self.PTRMS_cart, self.correct_C, self.TRM_1)
        print self.TRM_1
        print "original cartesian ptrm:", self.PTRMS_cart
        print "C (difference)", self.correct_C
        print "TRM_1", self.TRM_1, "rest of TRM", self.PTRMS_cart + self.correct_C
        for num, trm in enumerate(TRM_star):
            for n, v in enumerate(trm):
                self.assertAlmostEqual(self.TRM_star[num][n], v)

# so ref_TRM_star should be TRM_1 then PTRMS_cart + correct_C
# then find vector lengths to get x_star

    def test_delta_pal(self):
#        delta_pal, x_star, b_star = lib_ptrm.get_delta_pal(self.C, self.PTRMS_Cart, self.TRM_1) # probs also needs self.NRMS, or perhaps self.y_Arai
        pass
    

#    def test_delta_pal(self):
#        pass
#        ptrm - ptrm_checks = diffs
#        for num, diff in enumerate(diffs):
#        pass
        

if __name__ == "__main__":
    unittest.main()
