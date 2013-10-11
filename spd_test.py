#!/usr/bin/env python

import unittest
import numpy
import copy
import spd
import known_values


# could use this directly in here:
#import new_lj_thellier_gui_spd as tgs
#gui = tgs.Arai_GUI()
#thing = PintPars(gui.Data, '0238x6011044', 473., 623.)
#    specimens = spd.gui.Data.keys()
#    thing1 = spd.PintPars(spd.gui.Data, specimens[3], 523., 773.)
#    thing1.calculate_all_statistics()
#    thing2 = spd.PintPars(spd.gui.Data, specimens[4], 273., 798.)
#    thing2.calculate_all_statistics()
#    thing3 = spd.PintPars(spd.gui.Data, specimens[5], 598, 698)
#    thing3.calculate_all_statistics()



class CheckParams(unittest.TestCase):

    # no init
    print "starting CheckParams test"
    ref = "reference"
    obj = copy.deepcopy(spd.thing)
    obj_pars = obj.pars.copy()  # prevents it from changing, yay
    obj.calculate_all_statistics()
    obj_new_pars = obj.pars
    pre_calculation_pars = ['specimen_int_n', 'lab_dc_field']
    post_calculation_pars = ['vector_diffs_segment', 'delta_x_prime', 'partial_vds', 'B_anc', 'SCAT', 'specimen_int', 'specimen_fvds', 'specimen_b_beta', 'vector_diffs', 'specimen_YT', 'magic_method_codes', 'specimen_vds', 'specimen_int_n', 'centroid', 'max_diff', 'FRAC', 'GAP-MAX', 'y_prime', 'best_fit_circle', 'delta_y_prime', 'B_anc_sigma', 'B_lab', 'specimen_b_sigma', 'specimen_b', 'specimen_g', 'specimen_XT', 'specimen_f', 'specimen_k', 'specimen_q', 'lab_dc_field', 'specimen_w', 'x_prime', 'SSE', 'specimen_g_lim', 'R_2corr']  # remember to update this as you add stats

    def test_for_params_before(self):
        for par in self.pre_calculation_pars:
            self.assertIn(par, self.obj_pars.keys())
            
    def test_for_params_after(self):
#        """
#        check that calculate_all_statistics() generates all expected pars
#        """
        for par in self.post_calculation_pars:
            self.assertIn(par, self.obj_new_pars.keys())

    def test_for_extra_params(self):
        """
        check that calculate_all_statistics doesn't generate any unexpected pars
        """
        for par in self.obj_new_pars.keys():
            self.assertIn(par, self.post_calculation_pars)

    def test_params_are_not_empty(self):
        for par in self.obj_pars.keys():
            value = self.obj_pars[par]
           # print value
            self.assertIsNotNone(value)
            
class CheckInitialAttributeValues(unittest.TestCase):
    obj = copy.deepcopy(spd.thing)
    obj_attributes = {'s':obj.s, 'datablock': obj.datablock, 'x_Arai': obj.x_Arai, 'y_Arai': obj.y_Arai, 't_Arai': obj.t_Arai, 'x_Arai_segment': obj.x_Arai_segment, 'y_Arai_segment': obj.y_Arai_segment, "x_Arai_mean": obj.x_Arai_mean, "y_Arai_mean": obj.y_Arai_mean, "x_tail_check": obj.x_tail_check, 'y_tail_check': obj.y_tail_check, 'tail_checks_temperatures': obj.tail_checks_temperatures, 'tail_checks_starting_temperatures': obj.tail_checks_starting_temperatures, 'x_ptrm_check': obj.x_ptrm_check, 'y_ptrm_check': obj.y_ptrm_check, 'ptrm_checks_temperatures': obj.ptrm_checks_temperatures, 'ptrm_checks_starting_temperatures': obj.ptrm_checks_starting_temperatures, 'zijdblock': obj.zijdblock, 'z_temperatures': obj.z_temperatures, 'start': obj.start, 'end': obj.end, 'pars': obj.pars, 'specimen_Data': obj.specimen_Data, 'tmin': obj.tmin, 'tmax': obj.tmax, 'tmin_K': obj.tmin_K, 'tmax_K': obj.tmax_K} 
    known_values = known_values.values

    # test a bunch of the initial values against known expected values.  this indirectly tests get_data and such.  I don't know, maybe that will be enough.  Then I can do more thorough tests for the other stuff. 
    def test_name(self):
        self.assertEqual(self.obj.s, '0238x6011044')

    def test_known_values(self):
        for key, value in self.known_values.iteritems():  # goes through all values
            if type(value) == int or type(value) == float: # can't iterate over int type or float
               # print type(value)
                self.assertEqual(value, self.obj_attributes[key])
            elif value != None and type(value) != dict:
               # print type(value)
                for num, item in enumerate(value):
                    message = "%s: known value = %s; obj_attribute = %s" %(key, value[:150], self.obj_attributes[key][:150])
                    if type(item) == float or type(item) == numpy.float64:
                        self.assertEqual(round(self.obj_attributes[key][num], 8), round(item, 8), message)
                    else:
                        self.assertEqual(self.obj_attributes[key][num], item, message)


class CheckFrac(unittest.TestCase):

    print spd.thing.pars
    obj = copy.deepcopy(spd.thing)
    obj.pars['specimen_vds'] = 1
    obj.pars['vector_diffs_segment'] = [0.031295984827684566, 0.024151118628312387, 0.036482667142194579, 0.059128016249387697, 0.034082675643388412, 0.090581343759213978]
    
    def test_Frac_known_value(self):
        frac = self.obj.get_FRAC()
        print "frac", frac
        self.assertEqual(round(frac, 12), .27572180625)

    def test_Frac_with_negative_input(self): # this one is silly, but I was making a point.  so, whatever
        self.obj.pars['specimen_vds'] = 2
        self.obj.pars['vector_diffs_segment'] = [3., 5., -.1, 8.]
#        print "negative FRAC", self.obj.get_FRAC()
        self.assertRaises(ValueError, self.obj.get_FRAC)


class CheckVDSsequence(unittest.TestCase):

    obj = copy.deepcopy(spd.thing)
    obj.York_Regression()
    stuff = obj.get_vds()
#    stuff['vds'] = -.1


    def test_for_negative_values(self):
        for k, v in self.stuff.items():
            self.assertGreaterEqual(v, 0) # none of these stats can possibly be negative numbers


class CheckR_2corr(unittest.TestCase):
    
    obj = copy.deepcopy(spd.thing)
    R_2corr = obj.get_R_2corr()

    def testPositiveOutput(self):
        self.assertGreater(self.R_2corr, 0)

    def testSizeOutput(self): # not absolutely sure this is true.  but it seems like it has to be.  
        self.assertLess(self.R_2corr, 1)

class CheckR_det2(self):
    
    

    # make super simple test cases for these functions.  do it by hand.  i.e., y = [1, 2] y_prime = [1.5, 1.5]
    # also make test that it does raise an error when dividing by zero
    # put in checks to make sure you aren't dividing by zero (in several functions)

    # consider architecture.  possibility of a math library that is then imported into spd.  makes testing easier.  


if __name__ == "__main__":
    unittest.main()

