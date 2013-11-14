#!/usr/bin/env python                                                                                            

import numpy
import unittest
import copy
import spd
import lib_tail_check_statistics as lib_tail


#class CheckpTRMparams(unittest.TestCase): 

class CheckTailSelection(unittest.TestCase):

    tmax = 50
    y_Arai = [4, 3.5, 3, 2., 1., .5]
    t_Arai = [10, 20, 30, 40, 50, 60]
    y_tail = [3.5, 3., 2.5, 1.1, .6 ]
    tail_temps = [10, 30, 40, 50, 60]
    obj = copy.deepcopy(spd.thing)

    def test_n_tail(self):
        ref_n_tail = 4
#        ref_tail_section.  specifying this is not neacessary if we always start at the first tail check
        n_tail = lib_tail.get_n_tail(self.tmax, self.tail_temps)
        self.assertEqual(ref_n_tail, n_tail)

    def test_n_tail_real_data(self):
        ref_n_tail = 4
        n_tail = self.obj.get_n_tail()
        self.assertEqual(ref_n_tail, n_tail)
        
              
    
