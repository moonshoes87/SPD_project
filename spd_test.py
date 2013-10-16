#!/usr/bin/env python

import unittest
import numpy
import copy
import spd
import known_values
import lj_library as lib


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
    post_calculation_pars = ['vector_diffs_segment', 'delta_x_prime', 'partial_vds', 'B_anc', 'SCAT', 'specimen_int', 'specimen_fvds', 'specimen_b_beta', 'vector_diffs', 'specimen_YT', 'specimen_vds', 'specimen_int_n', 'centroid', 'max_diff', 'FRAC', 'GAP-MAX', 'y_prime', 'best_fit_circle', 'delta_y_prime', 'B_anc_sigma', 'B_lab', 'specimen_b_sigma', 'specimen_b', 'specimen_g', 'specimen_XT', 'specimen_f', 'specimen_k', 'specimen_q', 'lab_dc_field', 'specimen_w', 'x_prime', 'SSE', 'specimen_g_lim', 'R_corr2', 'R_det2', 'count_IZ', 'count_ZI', 'x_err', 'y_err', 'x_tag', 'y_tag', '_SCAT', 'Z']  # remember to update this as you add stats.  removed magic_method_codes...

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
    known_values = known_values.initial_values

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


class CheckYorkRegression(unittest.TestCase):
    
    obj = copy.deepcopy(spd.thing)
    obj.York_Regression()
    known_values = known_values.York_Regression_values
    obj_pars = obj.pars
    
    def test_York_Regression(self):
        for key, value in self.known_values.iteritems():  # goes through all values
            if type(value) == int or type(value) == float: # can't iterate over int type or float
               # print type(value)
                self.assertEqual(value, self.obj_pars[key])
            elif value != None and type(value) != dict:
               # print type(value)
                for num, item in enumerate(value):
                    message = "%s: known value = %s; obj_attribute = %s" %(key, value[:150], self.obj_pars[key][:150])
                    if type(item) == float or type(item) == numpy.float64:
                        self.assertEqual(round(self.obj_pars[key][num], 8), round(item, 8), message)
                    else:
                        self.assertEqual(self.obj_pars[key][num], item, message)
        



class CheckVDSsequence(unittest.TestCase): # NOT DONE come back to this

    obj = copy.deepcopy(spd.thing)
    obj.York_Regression()
    stuff = obj.get_vds()
#    stuff['vds'] = -.1


    def test_for_negative_values(self):
        for k, v in self.stuff.items():
            self.assertGreaterEqual(v, 0) # none of these stats can possibly be negative numbers



class CheckFrac(unittest.TestCase): # basically good to go

    #print spd.thing.pars
    obj = copy.deepcopy(spd.thing)
    obj.pars['specimen_vds'] = 2
    obj.pars['vector_diffs_segment'] = [1., 1.5, 3.]
    # trying it a little differently:
    vds = 2
    vector_diffs_segment = [1., 1.5, 3.]

    def test_FRAC(self):
        frac = lib.get_FRAC(self.vds, self.vector_diffs_segment)
        self.assertEqual(frac, 2.75)

    def test_FRAC_with_zero_vds(self):
        self.assertRaises(ValueError, lib.get_FRAC, 0, self.vector_diffs_segment)

    def test_FRAC_with_negative_input(self):
        self.assertRaises(ValueError, lib.get_FRAC, 1, [1., -.1, 2.])

    def test_lib_vs_actual(self): # checks that it is calculated the same in lib and in practice
        self.obj.pars['vds'] = 2
        self.obj.pars['vector_diffs_segment'] = [1., 1.5, 3.]
        frac = lib.get_FRAC(self.vds, self.vector_diffs_segment)
        obj_frac = self.obj.get_FRAC()
        self.assertEqual(frac, obj_frac)

    

class CheckR_corr2(unittest.TestCase):
    obj = copy.deepcopy(spd.thing)
    R_corr2 = obj.get_R_corr2()
    x_segment, y_segment = numpy.array([1., 5.]), numpy.array([0., 2.])
    x_avg = sum(x_segment) / len(x_segment)
    y_avg = sum(y_segment) / len(y_segment)

    def testPositiveOutput(self):
        """should produce positive output"""
        self.assertGreater(self.R_corr2, 0)

    def testSizeOutput(self): # not absolutely sure this is true.  but it seems like it has to be.  
        """output should be less than 1"""
        self.assertLess(self.R_corr2, 1)

    def testSimpleInput(self):
        """should produce expected output with simple input"""
        r = lib.get_R_corr2(self.x_avg, self.y_avg, self.x_segment, self.y_segment)
        self.assertEqual(.5, r)
        
    def testDivideByZero(self):
        """should raise ValueError when attempting to divide by zero"""
        self.assertRaises(ValueError, lib.get_R_corr2, 1., 1., numpy.array([1.]), numpy.array([1.]))

class CheckR_det2(unittest.TestCase): # acceptable working test
    y_segment = [1., 2.5, 3]
    y_avg = 2.
    y_prime = [1., 2., 3.]
    
    def test_simple_input(self):
        result = lib.get_R_det2(self.y_segment, self.y_avg, self.y_prime)
        self.assertEqual((1 - .25/2.25), result)
    


class CheckZigzag(unittest.TestCase):
    x = [1., 2., 3.]
    y = [0., 4., 5.]
    y_int = 5.
    x_int = 1.
    n = len(x)
    reference_b_wiggle = [5., .5, 0.]
    slope = 1.2
    Z = 8.8
    Z_star = 88.
    obj = copy.deepcopy(spd.thing)
    obj.x_Arai, obj.y_Arai = x, y
    obj.pars['specimen_YT'], obj.pars['specimen_XT'] = y_int, x_int
    obj.pars['specimen_b'], obj.n = slope, n
    
    def testWiggleB(self):
        for num, b in enumerate(self.reference_b_wiggle):
            result = lib.get_b_wiggle(self.x[num], self.y[num], self.y_int)
            self.assertAlmostEqual(result, b)

    def testZ(self):
        # x_segment, y_segment, x_int, y_int, slope
        result = lib.get_Z(self.x, self.y, self.x_int, self.y_int, self.slope)
        self.assertAlmostEqual(self.Z, result)

    def testPintParsZ(self):
        result = self.obj.get_Z()
        self.assertAlmostEqual(self.Z, result)

    def testZStar(self):
        result = lib.get_Zstar(self.x, self.y, self.x_int, self.y_int, self.slope, self.n)
        self.assertAlmostEqual(self.Z_star, result)

    def testPintParsZstar(self):
        result = self.obj.get_Zstar()
        self.assertAlmostEqual(self.Z_star, result)

class IZZI_MD(unittest.TestCase):
    points = numpy.array([1., 2., 3.])
    norm = 4.
    ref_normed_points = numpy.array([.25, .5, .75])
    x = numpy.array([4, 6, 12])
    y = numpy.array([8, 4, 2])
    norm_x = lib.get_normed_points(x, norm)
    norm_y = lib.get_normed_points(y, norm)

    L1 = numpy.sqrt(1.25)
    L2 = numpy.sqrt(2.5)
    L3 = numpy.sqrt(6.25)

    def testPointNorming(self): # satisfactory
        result = lib.get_normed_points(self.points, self.norm)
        for num, point in enumerate(result):
            self.assertAlmostEqual(self.ref_normed_points[num], point)

    def testTriangleSides(self): #
        result = lib.get_triangle_sides(self.norm_x, self.norm_y)
        line1, line2, line3 = result['L1'], result['L2'], result['L3']
        lines = [line1, line2, line3]
        ref_lines = [self.L1, self.L2, self.L3]
        for num, line in enumerate(lines):
            self.assertAlmostEqual(line, ref_lines[num])
#    get_triangle_coordinates  -- divide up points into appropriate ZI IZ ZI triangles
#    get_triangle_sides -- get the sides
#    get_triangle_area -- get the area
#    get_triangle_sign -- find if it is negative or positive


if __name__ == "__main__":
    unittest.main()

