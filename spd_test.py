#!/usr/bin/env python

import unittest
import spd
import known_values


#class ToRomanBadInput(unittest.TestCase):                            
#    def testTooLarge(self):                                          
#        """toRoman should fail with large input"""                   
#        self.assertRaises(roman.OutOfRangeError, roman.toRoman, 4000)

class CheckParams(unittest.TestCase):

    # no init
    print "starting CheckParams test"
    ref = "reference"
    obj = spd.thing
# could use this directly in here:
#import new_lj_thellier_gui_spd as tgs
#gui = tgs.Arai_GUI()
#thing = PintPars(gui.Data, '0238x6011044', 473., 623.)
    obj_pars = obj.pars.copy()  # prevents it from changing, yay
    obj.calculate_all_statistics()
    obj_new_pars = obj.pars
    pre_calculation_pars = ['specimen_int_n', 'lab_dc_field']
    post_calculation_pars = ['vector_diffs_segment', 'delta_x_prime', 'partial_vds', 'B_anc', 'SCAT', 'specimen_int', 'specimen_fvds', 'specimen_b_beta', 'vector_diffs', 'specimen_YT', 'magic_method_codes', 'specimen_vds', 'specimen_int_n', 'centroid', 'max_diff', 'FRAC', 'GAP-MAX', 'y_prime', 'best_fit_circle', 'delta_y_prime', 'B_anc_sigma', 'B_lab', 'specimen_b_sigma', 'specimen_b', 'specimen_g', 'specimen_XT', 'specimen_f', 'specimen_k', 'specimen_q', 'lab_dc_field', 'specimen_w', 'x_prime', 'SSE', 'specimen_g_lim']  # remember to update this as you add stats

    def test_for_params_before(self):
        for par in self.pre_calculation_pars:
            self.assertIn(par, self.obj_pars.keys())
            
    def test_for_params_after(self):
#        """
#        check that calculate_all_statistics() generates all expected pars
#        """
        for par in self.post_calculation_pars:
            self.assertIn(par, self.obj_new_pars.keys())

    def test_params_are_not_empty(self):
        for par in self.obj_pars.keys():
            value = self.obj_pars[par]
           # print value
            self.assertIsNotNone(value)
            
class CheckInitialValues(unittest.TestCase):
    obj = spd.thing
    stuff = [obj.s, obj.datablock, obj.x_Arai, obj.y_Arai, obj.t_Arai, obj.x_Arai_segment, obj.y_Arai_segment, obj.x_tail_check, obj.y_tail_check, obj.tail_checks_temperatures, obj.tail_checks_starting_temperatures, obj.x_ptrm_check, obj.y_ptrm_check, obj.ptrm_checks_temperatures, obj.ptrm_checks_starting_temperatures, obj.zijdblock, obj.z_temperatures, obj.start, obj.end, obj.pars, obj.specimen_Data, obj.tmin, obj.tmax, obj.tmin_K, obj.tmax_K] 
    known_values = known_values.values

    # test a bunch of the initial values against known expected values.  this indirectly tests get_data and such.  I don't know, maybe that will be enough.  Then I can do more thorough tests for the other stuff. 

    def test_name(self):
        self.assertEqual(self.obj.s, '0238x6011044')

    def test_known_values(self):
        n = 0
        for k, v in self.known_values.iteritems():  # not in order
            print k, v


# this is not syntactically correct, but you get the idea
#    def test_known_values(self):
#        for key, val in known_values:
#            self.assertEqual(self.obj.key == val)



#    def test_name_length(self):
#        print self.obj.pars
#        print unittest
#        self.assertEqual(10, 10)
#        self.assertEqual(len(self.obj.s), 12) 
        
#    def test_thing_is_pie(self):
#        self.assertRaises(ValueError, self.obj.York_Regression)

#    def test_SCAT(self):
#        print obj_pars
#        self.assertEqual(obj_pars['SCAT'], False)

#    def test_should_start_without_params(self):
#        self.assertFalse(1 == 2)

if __name__ == "__main__":
#    new_test = TestSpecimenValues(thing())
    unittest.main()

