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
# Author: Lori Jonestrask / Ron Shaar
#
# Initial revision: September 2013: 
#
#============================================================================================

import sys
#import pylab
#import scipy
from scipy import * 
#import os
import lib.lib_arai_plot_statistics as lib_arai
import lib.lib_curvature as lib_k
import lib.lib_directional_statistics as lib_direct
import lib.lib_ptrm_statistics as lib_ptrm
import lib.lib_tail_check_statistics as lib_tail
import lib.lib_additivity_check_statistics as lib_add


# Data{} is a dictionary sorted by specimen name
# i.e. Data[specine_a]={} ; Data[specine_b]={} etc.
# Each specimen data is sorted to th efollowing "blocks":

# Arai plot:
# Data[s]['x_Arai']=[] # a list of x_i 
# Data[s]['y_Arai']=[] # a list of y_i
# Data[s]['t_Arai']=[] # a list of temperatures (C? K?)
# Data[s]'steps_Arai']=[] a list of the "ZI/IZ" steps (IZ= blue ZI=red). i.e. ['ZI','IZ','ZI']

# NRMS  # temp K, dec, inc, moment, 0 or 1 (ZI=1, IZ=0)  # needs to convert to cartesian for vector operations
# PTRMS  # temp K, dec, inc, moment, 0 or 1 (ZI=1, IZ=0)

#PTRM_Checks # same same

# pTRM checks ("triangles")
#         Data[s]['x_ptrm_check']=[] # a list of x coordinates of pTRM checks
#         Data[s]['y_ptrm_check']=[] # a list of y coordinates of pTRM checks      
#         Data[s]['ptrm_checks_temperatures']=[] # a list of pTRM checks temperature 
#         Data[s]['x_ptrm_check_starting_point'] # a list of x coordinates of the point ehere the pTRM checks started from
#         Data[s]['y_ptrm_check_starting_point'] # a list of y coordinates of the point ehere the pTRM checks started from            
#         Data[s]['ptrm_checks_starting_temperatures']=[] # a list of temperatures from which the pTRM checks started from 


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
        #print "calling __init__ PintPars object"
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

        # AC checks
        self.x_add_check = self.specimen_Data['x_additivity_check']
        self.y_add_check = self.specimen_Data['y_additivity_check']
        self.add_checks_temperatures = self.specimen_Data['additivity_check_temperatures']
        self.add_checks_starting_temperatures = self.specimen_Data['additivity_check_starting_temperatures']
        self.x_add_check_starting_point = self.specimen_Data['x_additivity_check_starting_point']
        self.y_add_check_starting_point = self.specimen_Data['y_additivity_check_starting_point']
        
        # Data in Temp/Dec/Inc/Int format
        self.PTRMS = self.specimen_Data['PTRMS']
        self.NRMS = self.specimen_Data['NRMS']
        self.NRM = self.specimen_Data['NRM']
        self.PTRM_Checks = self.specimen_Data['PTRM_Checks']
        self.TAIL_Checks = self.specimen_Data['TAIL_Checks']
        self.AC_Diffs = self.specimen_Data['AC']
#        self.ADD_Checks = self.specimen_Data['ADD_Checks'] # removed this from new_lj_thellier_gui

        self.zijdblock=self.specimen_Data['zijdblock']        
        self.z_temperatures=self.specimen_Data['z_temp']

        self.start=self.t_Arai.index(tmin)
        self.end=self.t_Arai.index(tmax)

        self.pars={}

        self.pars['lab_dc_field']=self.specimen_Data['pars']['lab_dc_field']
        self.B_lab_dir = [self.specimen_Data['Thellier_dc_field_phi'], self.specimen_Data['Thellier_dc_field_theta'], 
                          self.specimen_Data['Thellier_dc_field_uT']]  # 

  #      self.pars['magic_method_codes']=Data[self.s]['pars']['magic_method_codes']
        self.pars['specimen_n']=self.end-self.start+1

 
        #LJ ADDING stats:
        self.steps_Arai = self.specimen_Data['steps_Arai']
        self.n = float(self.end-self.start+1)
        self.n_max = len(self.t_Arai)
        self.tmin = tmin 
        self.tmax = tmax
        self.tmin_K = tmin - 273. 
        self.tmax_K = tmax - 273
        self.x_Arai_segment = self.x_Arai[self.start:self.end+1] 
        self.y_Arai_segment = self.y_Arai[self.start:self.end+1]
        self.x_Arai_mean = mean(self.x_Arai_segment)
        self.y_Arai_mean = mean(self.y_Arai_segment)
        self.xy_Arai = lib_arai.get_xy_array(self.x_Arai, self.y_Arai)
        self.xy_Arai_segment = lib_arai.get_xy_array(self.x_Arai_segment, self.y_Arai_segment)


    def __repr__(self):
        return "PintPars object, specimen: {}, tmin_K: {}, tmax_K: {} ".format(self.s, self.tmin_K, self.tmax_K)

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
        return {'max_diff': data['max_diff'], 'vector_diffs': data['vector_diffs'], 'specimen_vds': data['specimen_vds'], 
                'specimen_fvds': data['specimen_fvds'], 'vector_diffs_segment': data['vector_diffs_segment'], 
                'partial_vds': data['partial_vds'], 'GAP-MAX': data['GAP-MAX']}

    def get_FRAC(self):
        vds = self.pars['specimen_vds']
        vector_diffs_segment = self.pars['vector_diffs_segment']
        FRAC = lib_arai.get_FRAC(vds, vector_diffs_segment)
        self.pars['FRAC'] = FRAC
        return FRAC

    def get_curve(self):
        x_Arai, y_Arai = self.x_Arai, self.y_Arai
        data = lib_k.AraiCurvature(x_Arai,y_Arai)
        self.pars['specimen_k'] = data[0]
        self.pars['SSE'] = data[3]
        return data[0], data[3]

    def get_SCAT(self):
        slope = self.pars['specimen_b'] 
        x_mean, y_mean = self.x_Arai_mean, self.y_Arai_mean
        x_Arai_segment, y_Arai_segment = self.x_Arai_segment, self.y_Arai_segment
        box = lib_arai.get_SCAT_box(slope, x_mean, y_mean)
    #    print "SCAT-box", box
        low_bound, high_bound, x_max, y_max = box[0], box[1], box[2], box[3]
        # getting SCAT points
        x_Arai_segment, y_Arai_segment = self.x_Arai_segment, self.y_Arai_segment
        tmin, tmax = self.tmin, self.tmax
        ptrm_checks_temps, ptrm_checks_starting_temps = self.ptrm_checks_temperatures, self.ptrm_checks_starting_temperatures
        x_ptrm_check, y_ptrm_check =  self.x_ptrm_check, self.y_ptrm_check
        tail_checks_temps, tail_checks_starting_temps =  self.tail_checks_temperatures, self.tail_checks_starting_temperatures
        x_tail_check, y_tail_check = self.x_tail_check, self.y_tail_check
        points = lib_arai.get_SCAT_points(x_Arai_segment, y_Arai_segment, tmin, tmax, 
                                          ptrm_checks_temps, ptrm_checks_starting_temps, 
                                          x_ptrm_check, y_ptrm_check, tail_checks_temps, 
                                          tail_checks_starting_temps, x_tail_check, y_tail_check)
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
        x_segment, y_segment = self.x_Arai_segment, self.y_Arai_segment
        x_int, y_int = self.pars['specimen_XT'], self.pars['specimen_YT']
        slope = self.pars['specimen_b']
        Z = lib_arai.get_Z(x_segment, y_segment, x_int, y_int, slope)
        self.pars['Z'] = Z
        return Z

    def get_Zstar(self):
        x_segment, y_segment = self.x_Arai_segment, self.y_Arai_segment
        x_int, y_int = self.pars['specimen_XT'], self.pars['specimen_YT']
        slope, n = self.pars['specimen_b'], self.n
        Zstar = lib_arai.get_Zstar(x_segment, y_segment, x_int, y_int, slope, n)
        self.pars['Zstar'] = Zstar
        return Zstar
                             
    def get_IZZI_MD(self):
        import lib.lib_IZZI_MD as lib_izzi
        if ('IZ' in self.steps_Arai):
            IZZI_MD = lib_izzi.get_IZZI_MD(self.x_Arai, self.y_Arai, self.steps_Arai, self.start, self.end)
            self.pars['IZZI_MD'] = IZZI_MD
            return IZZI_MD
        else:
            IZZI_MD = float('nan')
            self.pars['IZZI_MD'] = IZZI_MD
            return IZZI_MD

        
    # directional statistics begin here:

    def get_dec_and_inc(self):
        Dec_Anc, Inc_Anc, best_fit_Anc, tau_Anc, V_Anc, mass_center = lib_direct.get_dec_and_inc(self.zdata,
                self.t_Arai, self.tmin, self.tmax, anchored=True)
        Dec_Free, Inc_Free, best_fit_Free, tau_Free, V_Free, mass_center = lib_direct.get_dec_and_inc(self.zdata, 
                self.t_Arai, self.tmin, self.tmax, anchored=False)
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
        self.pars['MAD_Free'], self.pars['MAD_Anc'] = MAD_Free, MAD_Anc
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
        NRM_dev = lib_direct.get_NRM_dev(self.pars['DANG'], self.pars['zdata_mass_center'], self.pars['specimen_YT'])
        self.pars['NRM_dev'] = NRM_dev
        return NRM_dev

    def get_theta(self):
        b_lab_dir = [self.B_lab_dir[0], self.B_lab_dir[1], 1.]
        ChRM = self.pars['best_fit_vector_Free'] # WAS USING THIS
        ChRM = self.pars['best_fit_vector_Anc'] # THIS AGREES WITH GREIG's CODE
        theta = lib_direct.get_theta(b_lab_dir, ChRM)
        self.pars['theta'] = theta
        return theta

    def get_gamma(self):
        B_lab_dir = [self.B_lab_dir[0], self.B_lab_dir[1], 1.] # dir
       # ptrm_vector = [self.PTRMS[-1][1], self.PTRMS[-1][2], 1] # dir
        ind = self.t_Arai.index(self.tmax)
        alt_ptrm_dir = [self.PTRMS[ind][1], self.PTRMS[ind][2], self.PTRMS[ind][3] / self.specimen_Data['NRM']] # dir
        alt_ptrm_cart = lib_direct.dir2cart(alt_ptrm_dir)
        ptrm_dir = [self.PTRMS[-1][1], self.PTRMS[-1][2], self.PTRMS[-1][3] / self.specimen_Data['NRM']] # dir
        gamma = lib_direct.get_gamma(B_lab_dir, alt_ptrm_dir) # this agrees with Greig's code.
#        gamma = lib_direct.get_gamma(B_lab_dir, ptrm_dir) # what I have been using
        self.pars['gamma'] = gamma
        return gamma


# ptrm statistics begin here
    
    def get_n_ptrm(self):
        tmin, tmax = self.tmin, self.tmax
        ptrm_temps, ptrm_starting_temps = self.ptrm_checks_temperatures, self.ptrm_checks_starting_temperatures
        n, steps = lib_ptrm.get_n_ptrm(tmin, tmax, ptrm_temps, ptrm_starting_temps)
#        print "n_ptrm", n, "steps", steps
        self.pars['n_ptrm'] = n
        self.pars['ptrm_checks_included_temps'] = steps
        
    def get_max_ptrm_check(self):
        ptrm_checks_included_temps = self.pars['ptrm_checks_included_temps']
        ptrm_checks = self.ptrm_checks_temperatures
        ptrm_x = self.x_ptrm_check
        x_Arai, t_Arai = self.x_Arai, self.t_Arai
        ptrm_checks, max_ptrm_check, sum_ptrm_checks, check_percent, sum_abs_ptrm_checks = lib_ptrm.get_max_ptrm_check(ptrm_checks_included_temps, 
            ptrm_checks, ptrm_x, t_Arai, x_Arai)
        self.pars['ptrm_checks'] = ptrm_checks
        self.pars['max_ptrm_check_percent'] = check_percent
        self.pars['max_ptrm_check'] = max_ptrm_check
        self.pars['sum_ptrm_checks'] = sum_ptrm_checks
        self.pars['sum_abs_ptrm_checks'] = sum_abs_ptrm_checks
        return max_ptrm_check

    def get_delta_CK(self):
#        def get_delta_CK(max_ptrm_check, x_int):
        delta_CK = lib_ptrm.get_delta_CK(self.pars['max_ptrm_check'], self.pars['specimen_XT'])
        self.pars['delta_CK'] = delta_CK
        return delta_CK

    def get_DRAT(self):
#        def get_DRAT(delta_y_prime, delta_x_prime, max_ptrm_check):
        DRAT, L = lib_ptrm.get_DRAT(self.pars['delta_x_prime'], self.pars['delta_y_prime'], self.pars['max_ptrm_check'])
        self.pars['DRAT'] = DRAT
        self.pars['length_best_fit_line'] = L
        return DRAT

    def get_max_DEV(self):
        max_DEV = lib_ptrm.get_max_DEV(self.pars['delta_x_prime'], self.pars['max_ptrm_check'])
        self.pars['max_DEV'] = max_DEV
        return max_DEV

    def get_CDRAT(self):
        CDRAT, CDRAT_prime = lib_ptrm.get_CDRAT(self.pars['length_best_fit_line'], self.pars['sum_ptrm_checks'], 
                                                self.pars['sum_abs_ptrm_checks'])
        self.pars['CDRAT'], self.pars['CDRAT_prime'] = CDRAT, CDRAT_prime
        return CDRAT, CDRAT_prime

    def get_DRATS(self):
        DRATS, DRATS_prime = lib_ptrm.get_DRATS(self.pars['sum_ptrm_checks'], self.pars['sum_abs_ptrm_checks'], self.x_Arai, self.end)
        self.pars['DRATS'] = DRATS
        self.pars['DRATS_prime'] = DRATS_prime
        return DRATS

    def get_mean_DRAT(self):
        mean_DRAT, mean_DRAT_prime = lib_ptrm.get_mean_DRAT(self.pars['sum_ptrm_checks'], self.pars['sum_abs_ptrm_checks'], 
                                                            self.pars['n_ptrm'], self.pars['length_best_fit_line'])
        self.pars['mean_DRAT'] = mean_DRAT
        self.pars['mean_DRAT_prime'] = mean_DRAT_prime
        return mean_DRAT, mean_DRAT_prime

    def get_mean_DEV(self):
        mean_DEV, mean_DEV_prime = lib_ptrm.get_mean_DEV(self.pars['sum_ptrm_checks'], self.pars['sum_abs_ptrm_checks'], 
                                                         self.pars['n_ptrm'], self.pars['delta_x_prime'])
        self.pars['mean_DEV'] = mean_DEV
        self.pars['mean_DEV_prime'] = mean_DEV_prime

    def get_delta_pal(self): 
        ptrms_segment, checks_segment = lib_ptrm.get_segments(self.PTRMS, self.PTRM_Checks, self.tmax)
        delta_pal = lib_ptrm.get_full_delta_pal(self.PTRMS, self.PTRM_Checks, self.NRM, self.pars['y_err'], 
                                                self.y_Arai_mean, self.pars['specimen_b'], self.start, self.end)
        self.pars['delta_pal'] = delta_pal


        # tail check statistics

    def get_n_tail(self):
        if len(self.tail_checks_temperatures) > 0:
            n_tail = lib_tail.get_n_tail(self.tmax, self.tail_checks_temperatures)
        else:
            n_tail = 0
        self.pars['n_tail'] = n_tail
        return n_tail

    def get_max_tail_check(self):
        if len(self.y_tail_check) > 0:
            tail_check_max, tail_check_diffs = lib_tail.get_max_tail_check(self.y_Arai, self.y_tail_check, self.t_Arai, 
                                                                           self.tail_checks_temperatures, self.pars['n_tail'])
        else:
            tail_check_max, tail_check_diffs = 0, [0]
        self.pars['tail_check_max'], self.pars['tail_check_diffs'] = tail_check_max, tail_check_diffs
        return tail_check_max, tail_check_diffs

    def get_DRAT_tail(self):
        DRAT_tail = lib_tail.get_DRAT_tail(self.pars['tail_check_max'], self.pars['length_best_fit_line'])
        self.pars['DRAT_tail'] = DRAT_tail
        return DRAT_tail

    def get_delta_TR(self):
    #    def get_delta_TR(tail_check_max, y_int)
        delta_TR = lib_tail.get_delta_TR(self.pars['tail_check_max'], self.pars['specimen_YT'])
        self.pars['delta_TR'] = delta_TR
        return delta_TR

    def get_MD_VDS(self):
        MD_VDS = lib_tail.get_MD_VDS(self.pars['tail_check_max'], self.pars['specimen_vds'])
        self.pars['MD_VDS'] = MD_VDS
        return MD_VDS


    # additivity check statistics start here

    def get_n_add(self):
        n_add = lib_add.get_n_add(self.add_checks_temperatures, self.add_checks_starting_temperatures, self.tmax)[1]
        self.pars['n_add'] = n_add
        return n_add

    def get_delta_AC(self):
#        print "doing delta AC"
        delta_AC = lib_add.get_delta_AC(self.pars['n_add'], self.AC_Diffs, self.pars['specimen_XT'])
        self.pars['delta_AC'] = delta_AC
        return delta_AC
      
    def arai_plot_statistics(self):
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

    def directional_statistics(self):
        self.York_Regression()
        self.get_vds()
        #
        self.get_dec_and_inc()
        self.get_MAD()
        self.get_alpha()
        self.get_DANG()
        self.get_NRM_dev()
        self.get_theta() # not necessarily done
        self.get_gamma() # ditto        


    def ptrm_check_statistics(self):# ptrm check statistics
        self.York_Regression()
        self.get_vds()
        #
        self.get_n_ptrm()
        self.get_max_ptrm_check()
        self.get_delta_CK()
        self.get_DRAT()
        self.get_max_DEV()
        self.get_CDRAT()
        self.get_DRATS()
        self.get_mean_DRAT()
        self.get_mean_DEV()
        self.get_delta_pal()

    def tail_check_statistics(self):
        # tail check statistics
        self.York_Regression()
        self.get_vds()
        self.get_n_ptrm()
        self.get_max_ptrm_check()
        self.get_DRAT()
        #
        self.get_n_tail()
        self.get_max_tail_check()
        self.get_DRAT_tail()
        self.get_delta_TR()
        self.get_MD_VDS()

    def additivity_check_statistics(self):
        self.get_n_add()
        self.get_delta_AC()


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
        self.get_theta() 
        self.get_gamma()
        # ptrm check statistics
        self.get_n_ptrm()
        self.get_max_ptrm_check()
        self.get_delta_CK()
        self.get_DRAT()
        self.get_max_DEV()
        self.get_CDRAT()
        self.get_DRATS()
        self.get_mean_DRAT()
        self.get_mean_DEV()
        self.get_delta_pal()
        # tail check statistics
        self.get_n_tail()
        self.get_max_tail_check()
        self.get_DRAT_tail()
        self.get_delta_TR()
        self.get_MD_VDS()
        # additivity check statistics
        self.get_n_add()
        self.get_delta_AC()
        print "done with calculate_all_statistics"


# K temps: [0.0, 100.0, 150.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0]
# C temps: [273, 373.0, 423.0, 473.0, 498.0, 523.0, 548.0, 573.0, 598.0, 623.0, 648.0, 673.0, 698.0, 723.0, 748.0, 773.0, 798.0, 823.0]
import new_lj_thellier_gui_spd as tgs

gui = tgs.Arai_GUI()
#gui2 = tgs.Arai_GUI('new_magic_measurements.txt')
gui3 = tgs.Arai_GUI('consistency_tests/Bowles_etal_2006_magic_measurements.txt')
#specimens = gui2.Data.keys()
thing = PintPars(gui.Data, '0238x6011044', 473., 623.)
thing.calculate_all_statistics()

#mat_thing = PintPars(gui3.Data, 'AL2360-1b', 273., 498.)
#mat_thing.calculate_all_statistics()
#[  273, 423.0, 473.0, 523.0, 573.0, 623.0, 673.0]
#[   0.,  150.,  200.,  250.,  300.,  350.,  400.]

spec = PintPars(gui3.Data, 'AL3031-1d', 273., 498.)
spec2 = PintPars(gui3.Data, 'AL2360-1d', 273, 498.)
spec.calculate_all_statistics()
spec2.calculate_all_statistics()

#thing1 = PintPars(gui.Data, specimens[3], 523., 773.)
#thing = PintPars(gui.Data,  '0238x6011044', 273., 798.)

#thing = PintPars(gui.Data, specimens[3], 523., 773.)
#thing = PintPars(gui.Data, specimens[4], 273., 798.)
#thing = PintPars(gui.Data, specimens[2], 273., 773.)


