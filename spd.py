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



class PintPars(object):
    def __init__(self, Data,specimen_name,tmin,tmax):
        self.Data=Data
        self.s=specimen_name
        self.datablock = Data[self.s]['datablock']

        self.x_Arai=Data[self.s]['x_Arai']
        self.y_Arai=Data[self.s]['y_Arai']
        self.t_Arai=Data[self.s]['t_Arai']

        self.x_tail_check=Data[self.s]['x_tail_check']
        self.y_tail_check=Data[self.s]['y_tail_check']

        self.zijdblock=Data[self.s]['zijdblock']        
        self.z_temperatures=Data[self.s]['z_temp']

        self.start=self.t_Arai.index(tmin)
        self.end=self.t_Arai.index(tmax)

        self.pars={}
        self.pars['lab_dc_field']=Data[self.s]['pars']['lab_dc_field']
        self.pars['magic_method_codes']=Data[self.s]['pars']['magic_method_codes']

        self.pars['specimen_int_n']=self.end-self.start+1

    def York_Regression(self):
        #-------------------------------------------------
        # York regresssion (York, 1967) following Coe et al.(1978)
        # calculate statsitics: n,b,b_sigma_,f,fvds,specimen_int,g,q
        # modified from pmag.py
        #-------------------------------------------------               

        x_Arai_segment= self.x_Arai[self.start:self.end+1]
        y_Arai_segment= self.y_Arai[self.start:self.end+1]

        x_Arai_mean=mean(x_Arai_segment)
        y_Arai_mean=mean(y_Arai_segment)

        # equations (2),(3) in Coe (1978) for b, sigma
        n=self.end-self.start+1
        x_err=x_Arai_segment-x_Arai_mean
        y_err=y_Arai_segment-y_Arai_mean

        # York b
        york_b=-1* sqrt( sum(y_err**2) / sum(x_err**2) )

        # york sigma
        york_sigma= sqrt ( (2 * sum(y_err**2) - 2*york_b*sum(x_err*y_err)) / ( (n-2) * sum(x_err**2) ) )

        # beta  parameter                
        beta_Coe=abs(york_sigma/york_b)

        # y_T is the intercept of the extrepolated line
        # through the center of mass (see figure 7 in Coe (1978))
        y_T = y_Arai_mean - york_b* x_Arai_mean

        # calculate the extarplated data points for f and fvds
        # (see figure 7 in Coe (1978))

        x_tag=(y_Arai_segment - y_T ) / york_b
        y_tag=york_b*x_Arai_segment + y_T

        # intersect of the dashed square and the horizontal dahed line  next to delta-y-5 in figure 7, Coe (1978)
        x_prime=(x_Arai_segment+x_tag) / 2
        y_prime=(y_Arai_segment+y_tag) / 2

        f_Coe=abs((y_prime[0]-y_prime[-1])/y_T)

        f_vds=abs((y_prime[0]-y_prime[-1])/self.Data[self.s]['vds'])

        g_Coe= 1 - (sum((y_prime[:-1]-y_prime[1:])**2) / sum((y_prime[:-1]-y_prime[1:]))**2 )

        q_Coe=abs(york_b)*f_Coe*g_Coe/york_sigma


        count_IZ= self.Data[self.s]['steps_Arai'].count('IZ')
        count_ZI= self.Data[self.s]['steps_Arai'].count('ZI')
        if count_IZ >1 and count_ZI >1:
            self.pars['magic_method_codes']="LP-PI-BT-IZZI"
        elif count_IZ <1 and count_ZI >1:
            self.pars['magic_method_codes']="LP-PI-ZI"
        elif count_IZ >1 and count_ZI <1:
            self.pars['magic_method_codes']="LP-PI-IZ"            
        else:
            self.pars['magic_method_codes']=""
            
        self.pars["specimen_b"]=york_b
        self.pars["specimen_int"]=-1*self.pars['lab_dc_field']*self.pars["specimen_b"]
        self.pars["specimen_YT"]=y_T       
        self.pars["specimen_b_sigma"]=york_sigma
        self.pars["specimen_b_beta"]=beta_Coe
        self.pars["specimen_f"]=f_Coe
        self.pars["specimen_fvds"]=f_vds
        self.pars["specimen_g"]=g_Coe
        self.pars["specimen_q"]=q_Coe
        self.pars['magic_method_codes']+=":IE-TT"
        if 'x_ptrm_check' in self.Data[self.s].keys():
            if len(self.Data[self.s]['x_ptrm_check'])>0:
                self.pars['magic_method_codes']+=":LP-PI-ALT-PTRM"
        if 'x_tail_check' in self.Data[self.s].keys():
            if len(self.Data[self.s]['x_tail_check'])>0:
                self.pars['magic_method_codes']+=":LP-PI-BT-MD"        


    def calculate_all_statistics(self):
        self.York_Regression()
        print self.pars
        
        


