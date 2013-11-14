#!/usr/bin/env python                                                                                            

import numpy
import unittest
import copy
import spd
import lib_tail_check_statistics as lib_tail


#class CheckpTRMparams(unittest.TestCase): 

class CheckTailSelection(unittest.TestCase):

    tmax = 50
    y_Arai = [4, 3.9, 3.5, 3, 2., 1.]
    t_Arai = [10, 20, 30, 40, 50, 60]
    y_tail = [3.5, 4.2, 3., 1.4, .6 ]
    tail_temps = [10, 30, 40, 50, 60]
    obj = copy.deepcopy(spd.thing)
    ref_n_tail = 4

#                [4,    3.5, 3,  2.,   1., .5]
#-               [3.5,  4.2, 3., 1., .6 ]
    ref_tail_check_max = .7
    ref_diffs = [ .5, -.7,  0,  .6]

    def test_n_tail(self):
        ref_n_tail = 4
#        ref_tail_section.  specifying this is not neacessary if we always start at the first tail check.  need to check that this is so
        n_tail = lib_tail.get_n_tail(self.tmax, self.tail_temps)
        self.assertEqual(ref_n_tail, n_tail)

    def test_n_tail_real_data(self):
        n_tail = self.obj.get_n_tail()
        self.assertEqual(self.ref_n_tail, n_tail)
              
    def test_max_tail_check(self):
        tail_check_max, tail_check_diffs = lib_tail.get_max_tail_check(self.y_Arai, self.y_tail, self.t_Arai, self.tail_temps, self.ref_n_tail)
        for num, diff in enumerate(tail_check_diffs):
            self.assertAlmostEqual(self.ref_diffs[num], diff)
        self.assertAlmostEqual(self.ref_tail_check_max, tail_check_max)
        
                      
        
