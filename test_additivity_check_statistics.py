#!/usr/bin/env python

import unittest
import lib_additivity_check_statistics as lib_add


#class CheckTailSelection(unittest.TestCase):

class CheckAdditivity(unittest.TestCase):

    temps =          [1, 2, 3, 4, 5]
    starting_temps = [2, 3, 4, 5, 6]
    tmax = 4

    ref_incl_temps = [1, 2, 3]
    ref_n_add = 3

    x_add_check = [0, .2, .4, .5, .7]

    x_Arai =      [0, .1, .35, .5, .6,  .8,  1., 1.1]
    t_Arai =      [1,  2, 2.5,  3., 3.5, 4., 5., 6.]
    
    ptrm_star = ptrm[2, 1]: .1, ptrm[3, 2]: .4, ptrm[4,3]: .3, ptrm[5,4]: .2, ptrm[6,5]: .1
    # think I additionally need the check diffs.  what I have is just x check values.
    # I believe it will be: pTRM of start temperature - pTRM of additivity check step.  then I can compare those values.  
    # so we are comparing: ptrm(t=2) - ptrm(t=1) (second is heated from room temp), with ptrm(t=2) - ptrm(t=1) (second is demagnetized from first)



#    T_i < T_j
#    pTRM*(T_j, T_i) = pTRM(T_j, T_0) - pTRM(T_i, T_0)
    # means expected value of pTRM gained between i and j should be the same as the pTRM gained between 0 and j minus the pTRM gained between 0 and i
#AC_i,j = pTRM_star(T_j,T_i) − pTRM(T_j,T_i)    



    def test_data_selection(self):
        incl_temps, n_add = lib_add.get_n_add(self.temps, self.starting_temps, self.tmax)
        for num, temp in enumerate(incl_temps):
            print temp, self.ref_incl_temps[num]
            self.assertEqual(temp, self.ref_incl_temps[num])
        self.assertEqual(n_add, self.ref_n_add)


    def test_pTRM_star(self):
        pass

    def test_ACs(self):
        pass
# get AC check diffs
# think will be comparing x_Arai at that temp with x_Arai_add at that point

        
        

    def test_delta_AC(self):
        pass
