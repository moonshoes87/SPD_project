#!/usr/bin/env python

#============================================================================================
# LOG HEADER:
#============================================================================================
#
#
#  
#
#-------------------------
#
# Author: 
#
# Initial revision: September 2013: 
#
#============================================================================================

import sys,pylab,scipy,os
import lib_arai_plot_statistics as lib_arai
import lib_curvature as lib_k
import lib_directional_statistics as lib_direct
import lib_ptrm_statistics as lib_ptrm
from scipy import * 


# Data{} is a dictionary sorted by specimen name
# i.e. Data[specine_a]={} ; Data[specine_b]={} etc.
# Each specimen data is sorted to th efollowing "blocks":

# Arai plot:
# Data[s]['x_Arai']=[] # a list of x_i 
# Data[s]['y_Arai']=[] # a list of y_i
# Data[s]['t_Arai']=[] # a list of temperatures (C? K?)
# Data[s]'steps_Arai']=[] a list of the "ZI/IZ" steps (IZ= blue ZI=red). i.e. ['ZI','IZ','ZI']

# pTRM checks ("triangles")
#         Data[s]['x_ptrm_check']=[] # a list of x coordinates of pTRM checks
#         Data[s]['y_ptrm_check']=[] # a list of y coordinates of pTRM checks      
#         Data[s]['ptrm_checks_temperatures']=[] # a list of pTRM checks temperature 
#         Data[s]['x_ptrm_check_starting_point'] # a list of x coordinates of the point ehere the pTRM checks started from
#         Data[s]['y_ptrm_check_starting_point'] # a list of y coordinates of the point ehere the pTRM checks started from            #         Data[s]['ptrm_checks_starting_temperatures']=[] # a list of temperatures from which the pTRM checks started from 

# pTRM tail checks ("squares")
#        Data[s]['x_tail_check']
#        Data[s]['y_tail_check']
#        Data[s]['tail_check_temperatures']
#        Data[s]['x_tail_check_starting_point']
#        Data[s]['y_tail_check_starting_point']
#        Data[s]['tail_checks_starting_temperatures']


# Zijderveld plot

# Data[s]['zijdblock']=[ [], [], []] # a list of:
#                   [treatment,declination,inclination,intensity,ZI,measurement_flag,magic_instrument_codes]
#                    treatment: temperature steps
#
# Data[s]['zij_rotated']=[[x,y,z], [x,y,z] ,[x,y,z] ]: a list of the x,y,z coordinates of the rotated Zijderveld plot
#                                                      the rotated zijderveld plot is what is plotted in Thellier Gui.
# Data[s]['zdata'] -- same as zij_rotated, but not rotated


class PintPars(object):
    def __init__(self, Data,specimen_name,tmin,tmax):
        print "calling __init__ PintPars object"
        self.s=specimen_name
        self.specimen_Data=Data[self.s]
        self.datablock = self.specimen_Data['datablock']

        self.x_Arai=self.specimen_Data['x_Arai']
        self.y_Arai=self.specimen_Data['y_Arai']
        self.t_Arai=self.specimen_Data['t_Arai']

        self.zdata = self.specimen_Data['zdata'] # LJ add

        self.x_tail_check=self.specimen_Data['x_tail_check']
        self.y_tail_check=self.specimen_Data['y_tail_check']
        self.tail_checks_temperatures = self.specimen_Data['tail_check_temperatures']
#        self.x_tail_check_starting_point = self.specimen_Data['x_tail_check_starting_point']
#        self.y_tail_check_starting_point = self.specimen_Data['y_tail_check_starting_point']
        self.tail_checks_starting_temperatures = self.specimen_Data['tail_checks_starting_temperatures']

# pTRM checks ("triangles")
        self.x_ptrm_check = self.specimen_Data['x_ptrm_check'] # a list of x coordinates of pTRM checks
        self.y_ptrm_check = self.specimen_Data['y_ptrm_check'] # a list of y coordinates of pTRM checks      
        self.ptrm_checks_temperatures= self.specimen_Data['ptrm_checks_temperatures'] # a list of pTRM checks temperature 
#        Data[s]['x_ptrm_check_starting_point'] # a list of x coordinates of the point ehere the pTRM checks started from
#         Data[s]['y_ptrm_check_starting_point'] # a list of y coordinates of the point ehere the pTRM checks started from     
        self.ptrm_checks_starting_temperatures = self.specimen_Data['ptrm_checks_starting_temperatures'] # a list of temperatures from which the pTRM checks started from 

        self.PTRMS = self.specimen_Data['PTRMS']

        self.zijdblock=self.specimen_Data['zijdblock']        
        self.z_temperatures=self.specimen_Data['z_temp']

        # tmax, tmin, start, and end -- inclusive, or exclusive???
        self.start=self.t_Arai.index(tmin)
        self.end=self.t_Arai.index(tmax)

        self.pars={}

        self.pars['lab_dc_field']=self.specimen_Data['pars']['lab_dc_field']
        self.B_lab_vector = [self.specimen_Data['Thellier_dc_field_phi'], self.specimen_Data['Thellier_dc_field_theta'], self.specimen_Data['Thellier_dc_field_uT']]  # although generally last one can be abstracted as 1

  #      self.pars['magic_method_codes']=Data[self.s]['pars']['magic_method_codes']
        self.pars['specimen_int_n']=self.end-self.start+1
        self.stuff = ["s", "datablock", "x_Arai", "y_Arai", "t_Arai", "x_Arai_segment", "y_Arai_segment", "x_Arai_mean", "y_Arai_mean", "x_tail_check", "y_tail_check", "tail_checks_temperatures", "tail_checks_starting_temperatures", "x_ptrm_check", "y_ptrm_check", "ptrm_checks_temperatures", "ptrm_checks_starting_temperatures", "zijdblock", "z_temperatures", "start", "end", "pars", "specimen_Data", "tmin", "tmax", "tmin_K", "tmax_K", "steps_Arai", "xy_Arai", "xy_Arai_segment", "B_lab_vector", "PTRMS"] # needs to be updated
 
        #LJ ADDING stats:
        self.steps_Arai = self.specimen_Data['steps_Arai']
        self.n = float(self.end-self.start+1)  # n is the number of points involved
        self.n_max = len(self.t_Arai)  # gets total number of temperatures taken.  (p. 4, top)
        self.tmin = tmin # self-explanatory
        self.tmax = tmax
        self.tmin_K = tmin - 273. #273.15
        self.tmax_K = tmax - 273 #273.15
        self.x_Arai_segment = self.x_Arai[self.start:self.end+1]  # returns array of relevant x points
        self.y_Arai_segment = self.y_Arai[self.start:self.end+1]
        self.x_Arai_mean = mean(self.x_Arai_segment) # uses scipy mean function to get the mean of the x points
        self.y_Arai_mean = mean(self.y_Arai_segment)
        self.xy_Arai = lib_arai.get_xy_array(self.x_Arai, self.y_Arai)
        self.xy_Arai_segment = lib_arai.get_xy_array(self.x_Arai_segment, self.y_Arai_segment)


    def get_segments_and_means(self):
        pass # consider making this a real deal thing.  
        

    def York_Regression(self):
        x_segment, y_segment = self.x_Arai_segment, self.y_Arai_segment
        x_mean, y_mean = self.x_Arai_mean, self.y_Arai_mean
        n = self.n
        lab_dc_field = float(self.specimen_Data['lab_dc_field'])
        steps_Arai = self.specimen_Data['steps_Arai']
        data = lib_arai.York_Regression(x_segment, y_segment, x_mean, y_mean, n, lab_dc_field, steps_Arai)
        self.pars['x_err'] = data['x_err']
        self.pars['y_err'] = data['y_err']
        self.pars['x_tag'] = data['x_tag']
        self.pars['y_tag'] = data['y_tag']
        self.pars['specimen_b'] = data['specimen_b']
        self.pars['specimen_b_sigma'] = data['specimen_b_sigma']
        self.pars['specimen_b_beta'] = data['specimen_b_beta']
        self.pars['specimen_YT'] = data['y_int']
        self.pars['specimen_XT'] = data['x_int']
        self.pars['x_prime'] = data['x_prime']
        self.pars['y_prime'] = data['y_prime']
        self.pars['delta_x_prime'] = data['delta_x_prime']
        self.pars['delta_y_prime'] = data['delta_y_prime']
        self.pars['specimen_f'] = data['specimen_f']
        self.pars['specimen_g'] = data['specimen_g']
        self.pars['specimen_g_lim'] = data['specimen_g_lim']
        self.pars['specimen_q'] = data['specimen_q']
        self.pars['specimen_w'] = data['specimen_w']
        self.pars['count_IZ'] = data['count_IZ']
        self.pars['count_ZI'] = data['count_ZI']
        self.pars['B_lab'] = data['B_lab']  # think I don't need this, actually
        self.pars['B_anc'] = data['B_anc']
        self.pars['B_anc_sigma'] = data['B_anc_sigma']
        self.pars['specimen_int'] = data['specimen_int']
        return data

    # eventually add a function to change tmax and/or tmin.  must also change start, end, 

    def get_vds(self):
        zdata = self.zdata
        delta_y_prime = self.pars['delta_y_prime']
        start, end = self.start, self.end
        data = lib_arai.get_vds(zdata, delta_y_prime, start, end)
        self.pars['max_diff'] = data['max_diff']
        self.pars['vector_diffs'] = data['vector_diffs']
        self.pars['specimen_vds'] = data['specimen_vds']
        self.pars['specimen_fvds']= data['specimen_fvds']
        self.pars['vector_diffs_segment'] = data['vector_diffs_segment']
        self.pars['partial_vds'] = data['partial_vds']
        self.pars['GAP-MAX'] = data['GAP-MAX']
        return {'max_diff': data['max_diff'], 'vector_diffs': data['vector_diffs'], 'specimen_vds': data['specimen_vds'], 'specimen_fvds': data['specimen_fvds'], 'vector_diffs_segment': data['vector_diffs_segment'], 'partial_vds': data['partial_vds'], 'GAP-MAX': data['GAP-MAX']}


    def get_FRAC(self):
        vds = self.pars['specimen_vds']
        vector_diffs_segment = self.pars['vector_diffs_segment']
        FRAC = lib_arai.get_FRAC(vds, vector_diffs_segment)
        self.pars['FRAC'] = FRAC
        return FRAC


    def get_curve(self):
        x_Arai, y_Arai = self.x_Arai, self.y_Arai
        data = lib_k.AraiCurvature(x_Arai,y_Arai)
        # k, a, b, SSE
        self.pars['specimen_k'] = data[0]
        self.pars['SSE'] = data[3]
        return data[0], data[3]
#        data = lib_arai.get_curve(x_Arai, y_Arai)
#        self.pars['centroid'] = data['centroid']
#        self.pars['specimen_k'] = data['k']
#        self.pars['best_fit_circle'] = data['best_fit_circle']
#        self.pars['SSE'] = data['SSE']


    def get_SCAT(self):
        slope = self.pars['specimen_b'] #, slope_err, slope_beta = self.pars['specimen_b'], self.pars['specimen_b_sigma'], self.pars['specimen_b_beta']
#        x_int, y_int = self.pars['specimen_XT'], self.pars['specimen_YT']
#        beta_threshold = .1                                                                                  
        x_mean, y_mean = self.x_Arai_mean, self.y_Arai_mean
        x_Arai_segment, y_Arai_segment = self.x_Arai_segment, self.y_Arai_segment
        box = lib_arai.get_SCAT_box(slope, x_mean, y_mean)
#def get_SCAT_box(slope,  x_mean, y_mean, beta_threshold = .1):
    #    print "SCAT-box", box
        low_bound, high_bound, x_max, y_max = box[0], box[1], box[2], box[3]
        # getting SCAT points
        x_Arai_segment, y_Arai_segment = self.x_Arai_segment, self.y_Arai_segment
        tmin, tmax = self.tmin, self.tmax
        ptrm_checks_temps, ptrm_checks_starting_temps, x_ptrm_check, y_ptrm_check = self.ptrm_checks_temperatures, self.ptrm_checks_starting_temperatures, self.x_ptrm_check, self.y_ptrm_check
        tail_checks_temps, tail_checks_starting_temps, x_tail_check, y_tail_check = self.tail_checks_temperatures, self.tail_checks_starting_temperatures, self.x_tail_check, self.y_tail_check
        points = lib_arai.get_SCAT_points(x_Arai_segment, y_Arai_segment, tmin, tmax, ptrm_checks_temps, ptrm_checks_starting_temps, x_ptrm_check, y_ptrm_check, tail_checks_temps, tail_checks_starting_temps, x_tail_check, y_tail_check)
        # checking each point
        SCAT = lib_arai.get_SCAT(points, low_bound, high_bound, x_max, y_max)
        self.pars['SCAT'] = SCAT
        return SCAT
        

    def get_R_corr2(self):
        x_avg = self.x_Arai_mean
        y_avg = self.y_Arai_mean
        x_segment =self.x_Arai_segment
        y_segment =self.y_Arai_segment
        R_corr2 = lib_arai.get_R_corr2(x_avg, y_avg, x_segment, y_segment)
        self.pars['R_corr2'] = R_corr2
        return R_corr2


    def get_R_det2(self):
        y_segment = self.y_Arai_segment
        y_avg = self.y_Arai_mean
        y_prime = self.pars['y_prime']
        R_det2 = lib_arai.get_R_det2(y_segment, y_avg, y_prime)
        self.pars['R_det2'] = R_det2

    def get_Z(self):
        x_segment, y_segment = self.x_Arai, self.y_Arai
        x_int, y_int = self.pars['specimen_XT'], self.pars['specimen_YT']
        slope = self.pars['specimen_b']
        Z = lib_arai.get_Z(x_segment, y_segment, x_int, y_int, slope)
        self.pars['Z'] = Z
        return Z


    def get_Zstar(self):
        x_segment, y_segment = self.x_Arai, self.y_Arai
        x_int, y_int = self.pars['specimen_XT'], self.pars['specimen_YT']
        slope, n = self.pars['specimen_b'], self.n
        Zstar = lib_arai.get_Zstar(x_segment, y_segment, x_int, y_int, slope, n)
        self.pars['Zstar'] = Zstar
        return Zstar
                             


    def get_IZZI_MD(self): # not sure if this is right.  
#def get_IZZI_MD(x_Arai=x_arai, y_Arai=y_arai, steps_Arai=steps_arai):
        x_Arai = self.x_Arai
        y_Arai = self.y_Arai
        steps_Arai = self.steps_Arai
        IZZI_MD = lib_arai.get_IZZI_MD(x_Arai, y_Arai, steps_Arai)
        self.pars['IZZI_MD'] = IZZI_MD
        return IZZI_MD

        
    # directional statistics begin here:

    def get_dec_and_inc(self):
        print self.s
        print "-"
        Dec_Anc, Inc_Anc, best_fit_Anc, tau_Anc, V_Anc, mass_center = lib_direct.get_dec_and_inc(self.zdata, self.t_Arai, self.tmin, self.tmax, anchored=True)
        Dec_Free, Inc_Free, best_fit_Free, tau_Free, V_Free, mass_center = lib_direct.get_dec_and_inc(self.zdata, self.t_Arai, self.tmin, self.tmax, anchored=False)
        self.pars['Dec_Anc'], self.pars['Dec_Free'] = Dec_Anc, Dec_Free
        self.pars['Inc_Anc'], self.pars['Inc_Free'] = Inc_Anc, Inc_Free
        self.pars['best_fit_vector_Anc'], self.pars['best_fit_vector_Free'] = best_fit_Anc, best_fit_Free
        self.pars['tau_Anc'], self.pars['tau_Free'] = tau_Anc, tau_Free
        self.pars['V_Anc'], self.pars['V_Free'] = V_Anc, V_Free
        self.pars['zdata_mass_center'] = mass_center

        
    def get_MAD(self):
        MAD_Free = lib_direct.get_MAD(self.pars['tau_Free'])
        MAD_Anc = lib_direct.get_MAD(self.pars['tau_Anc'])
#        MAD_Lisa_Free = lib_direct.Lisa_get_MAD(self.pars['tau_Free'])
        self.pars['MAD_free'], self.pars['MAD_anc'] = MAD_Free, MAD_Anc
#        self.pars['Lisa_MAD'] = MAD_Anc_lisa
#        self.pars['MAD_lisa_free'] = MAD_Lisa_Free
        return {'MAD_Free': MAD_Free, 'MAD_Anc': MAD_Anc }
    
       
    def get_alpha(self): # need Int_Free and Int_Anc
        free = self.pars['best_fit_vector_Free']
        anc = self.pars['best_fit_vector_Anc']
        alpha = lib_direct.get_alpha(anc, free)
        self.pars['alpha'] = alpha

    def get_DANG(self):
        free = self.pars['best_fit_vector_Free']
        cm = self.pars['zdata_mass_center']
        DANG = lib_direct.get_angle_difference(free, cm)
        self.pars['DANG'] = DANG

    def get_NRM_dev(self):
#        get_NRM_dev(self.dang, self.X_avg, self.y_int)
        NRM_dev = lib_direct.get_NRM_dev(self.pars['DANG'], self.pars['zdata_mass_center'], self.pars['specimen_YT'])
        self.pars['NRM_dev'] = NRM_dev
        return NRM_dev

    def get_theta(self):
        b_lab_dir = [self.B_lab_vector[0], self.B_lab_vector[1], 1.]
        print self.s
        ChRM = self.pars['best_fit_vector_Free']
        print "in spd: b_lab_dir", b_lab_dir, "ChRM", ChRM
        theta = lib_direct.get_theta(b_lab_dir, ChRM)
        self.pars['theta'] = theta
        return theta

    def get_gamma(self):
        lab_vector = [self.B_lab_vector[0], self.B_lab_vector[1], 1.] # dir
        ptrm_vector = [self.PTRMS[-1][1], self.PTRMS[-1][2], 1] # dir
        gamma = lib_direct.get_gamma(lab_vector, ptrm_vector)
        self.pars['gamma'] = gamma
        return gamma

# ptrm statistics begin here
    
    def get_n_ptrm(self):
        tmin, tmax = self.tmin, self.tmax
        ptrm_temps, ptrm_starting_temps = self.ptrm_checks_temperatures, self.ptrm_checks_starting_temperatures
        n, steps = lib_ptrm.get_n_ptrm(tmin, tmax, ptrm_temps, ptrm_starting_temps)
#        print "n_ptrm", n, "steps", steps
        self.pars['n_ptrm'] = n
        self.pars['ptrm_checks_segment'] = steps
        
    def get_max_ptrm_check(self):
        ptrm_checks_segment = self.pars['ptrm_checks_segment']
        ptrm_checks = self.ptrm_checks_temperatures
        ptrm_x = self.x_ptrm_check
        x_Arai, t_Arai = self.x_Arai, self.t_Arai
        max_ptrm_check, sum_ptrm_checks, check_percent = lib_ptrm.get_max_ptrm_check(ptrm_checks_segment, ptrm_checks, ptrm_x, t_Arai, x_Arai)
        self.pars['max_ptrm_check_percent'] = check_percent
        self.pars['max_ptrm_check'] = max_ptrm_check
        self.pars['sum_ptrm_checks'] = sum_ptrm_checks
        return max_ptrm_check

    def get_delta_CK(self):
#        def get_delta_CK(max_ptrm_check, x_int):
        delta_CK = lib_ptrm.get_delta_CK(self.pars['max_ptrm_check'], self.pars['specimen_XT'])
        self.pars['delta_CK'] = delta_CK
        return delta_CK

    def get_DRAT(self):
#        def get_DRAT(delta_y_prime, delta_x_prime, max_ptrm_check):
        DRAT = lib_ptrm.get_DRAT(self.pars['delta_y_prime'], self.pars['delta_x_prime'], self.pars['max_ptrm_check'])
        self.pars['DRAT'] = DRAT
        return DRAT
                                 

    def calculate_all_statistics(self):
        print "calling calculate_all_statistics in spd.py"
        self.York_Regression()
        self.get_vds()
        self.get_FRAC()
        self.get_curve()
        self.get_SCAT()
        self.get_R_corr2()
        self.get_R_det2()
        self.get_Z()
        self.get_Zstar()
        self.get_IZZI_MD()
        # directional statistics
        self.get_dec_and_inc()
        self.get_MAD()
        self.get_alpha()
        self.get_DANG()
        self.get_NRM_dev()
        self.get_theta() # not necessarily done
        self.get_gamma() # ditto
        # ptrm check statistics
        self.get_n_ptrm()
        self.get_max_ptrm_check()
        self.get_delta_CK()
        self.get_DRAT()
        print "done with calculate_all_statistics"


# K temps: [0.0, 100.0, 150.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0]
# C temps: [273, 373.0, 423.0, 473.0, 498.0, 523.0, 548.0, 573.0, 598.0, 623.0, 648.0, 673.0, 698.0, 723.0, 748.0, 773.0, 798.0, 823.0]
import new_lj_thellier_gui_spd as tgs
gui = tgs.Arai_GUI()
thing = PintPars(gui.Data, '0238x6011044', 473., 623.)
specimens = gui.Data.keys()
thing1 = PintPars(gui.Data, specimens[3], 523., 773.)
#thing = PintPars(gui.Data,  '0238x6011044', 273., 798.)

#thing = PintPars(gui.Data, specimens[3], 523., 773.)
#thing = PintPars(gui.Data, specimens[4], 273., 798.)
#thing = PintPars(gui.Data, specimens[2], 273., 773.)
thing.calculate_all_statistics()

if False:
    gui = tgs.Arai_GUI()
    thing = PintPars(gui.Data, '0238x6011044', 473., 623.) 
    gui = tgs.Arai_GUI()
    specimens = gui.Data.keys()
    thing = PintPars(gui.Data, '0238x6011044', 473., 623.) 
    thing.calculate_all_statistics()
    thing1 = PintPars(gui.Data, specimens[3], 523., 773.)
    thing1.calculate_all_statistics()
    thing2 = PintPars(gui.Data, specimens[4], 273., 798.)
    thing2.calculate_all_statistics()
    thing3 = PintPars(gui.Data, specimens[5], 598, 698)
    thing3.calculate_all_statistics()
    thing4 = PintPars(gui.Data, specimens[2], 273., 773.)
    thing4.calculate_all_statistics()
    print "---"

