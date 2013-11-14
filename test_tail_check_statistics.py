#!/usr/bin/env python                                                                                            

import numpy
import unittest
import spd
import lib_tail_check_statistics as lib_tail


#class CheckpTRMparams(unittest.TestCase): 

class CheckTailSelection(unittest.TestCase):

    tmax = 50
    y_Arai = [4, 3.5, 3, 2., 1., .5]
    t_Arai = [10, 20, 30, 40, 50, 60]
    y_tail = [3.5, 3., 2.5, 1.1, .6 ]
    tail_temps = [10, 30, 40, 50, 60]
              
    
