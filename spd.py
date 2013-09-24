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
        print "calling __init__ PintPars object"
        self.s=specimen_name
        self.specimen_Data=Data[self.s]
        print self.s
        self.datablock = self.specimen_Data['datablock']

        self.x_Arai=self.specimen_Data['x_Arai']
        self.y_Arai=self.specimen_Data['y_Arai']
        self.t_Arai=self.specimen_Data['t_Arai']
#        print self.t_Arai 

        self.x_tail_check=self.specimen_Data['x_tail_check']
        self.y_tail_check=self.specimen_Data['y_tail_check']

        self.zijdblock=self.specimen_Data['zijdblock']        
        self.z_temperatures=self.specimen_Data['z_temp']

        # tmax, tmin, start, and end -- inclusive, or exclusive???
        self.start=self.t_Arai.index(tmin)
        self.end=self.t_Arai.index(tmax)

        print "tmin, tmax"
        print tmin, tmax
        print "start", "end"
        print self.start, self.end
        # name of object at end is p

        self.pars={}

        print "pars"
        print self.specimen_Data['pars']
        self.pars['lab_dc_field']=self.specimen_Data['pars']['lab_dc_field']
  #      self.pars['magic_method_codes']=Data[self.s]['pars']['magic_method_codes']
        # for some reason missing any magic_method_codes.  possibly these would have been incorporated into the data from rmag_anisotropy or something
        # magic_method codes are locked up in datablock, not actually extracted.  not sure if this happens somewhere else in thellier_gui or not
        # also, fix the weirdness of having to set the precise number for tmin and tmax
        self.pars['specimen_int_n']=self.end-self.start+1
        self.stuff = ["s", "datablock", "x_Arai", "y_Arai", "t_Arai", "x_tail_check", "y_tail_check", "zijdblock", "z_temperatures", "start", "end", "pars", "specimen_Data"] # needs to be updated
 
        #LJ ADDING stats:
        self.n_max = len(self.t_Arai)  # gets total number of temperatures taken.  (p. 4, top)
        self.tmin = tmin # self-explanatory
        self.tmax = tmax

    # eventually add a function to change tmax and/or tmin.  must also change start, end, 


    def York_Regression(self):
        #-------------------------------------------------
        # York regresssion (York, 1967) following Coe et al.(1978)
        # calculate statsitics: n,b,b_sigma_,f,fvds,specimen_int,g,q
        # modified from pmag.py
        #-------------------------------------------------               

        x_Arai_segment= self.x_Arai[self.start:self.end+1]  # returns array of relevant x points
        y_Arai_segment= self.y_Arai[self.start:self.end+1]

        x_Arai_mean=mean(x_Arai_segment) # uses scipy mean function to get the mean of the x points
        y_Arai_mean=mean(y_Arai_segment)

        # equations (2),(3) in Coe (1978) for b, sigma
        n=self.end-self.start+1  # n is the number of points involved
        x_err=x_Arai_segment-x_Arai_mean  # seems to subtract the mean from each number in the x array (???).  it is a scipy array
        y_err=y_Arai_segment-y_Arai_mean

        # York b
        york_b=-1* sqrt( sum(y_err**2) / sum(x_err**2) ) # averaged slope
        
        # could we not use self.n_max instead of n here?
        # york sigma
        york_sigma= sqrt ( (2 * sum(y_err**2) - 2*york_b* sum(x_err*y_err)) / ( (n-2) * sum(x_err**2) ) )

        # beta  parameter                
        beta_Coe=abs(york_sigma/york_b)  # absolute value of york sigma/york_b

        # y_T is the intercept of the extrepolated line
        # through the center of mass (see figure 7 in Coe (1978))
        y_T = y_Arai_mean - york_b* x_Arai_mean

        #LJ x_T is x intercept
        x_T = (-1 * y_T) / york_b # LJ added


        # calculate the extarplated data points for f and fvds
        # (see figure 7 in Coe (1978))

        x_tag=(y_Arai_segment - y_T ) / york_b # returns array of y points minus the y intercept, divided by slope
        y_tag=york_b*x_Arai_segment + y_T

#        self.pars['x_tag'] = x_tag # LJ add

        # intersect of the dashed square and the horizontal dahed line  next to delta-y-5 in figure 7, Coe (1978)
        x_prime=(x_Arai_segment+x_tag) / 2
        y_prime=(y_Arai_segment+y_tag) / 2
        self.pars['x_prime'] = x_prime #LJ add
        self.pars['y_prime'] = y_prime # LJ add

        delta_x_prime = abs(x_prime[-1] - x_prime[0]) #Lj add.  this is the TRM length of the best fit line
        delta_y_prime = abs(y_prime[-1] - y_prime[0]) # LJ add.  this is the NRM length of the best fit line
        self.pars['delta_x_prime'] = delta_x_prime
        self.pars['delta_y_prime'] = delta_y_prime

        f_Coe=abs((y_prime[0]-y_prime[-1])/y_T)  # same as 'f' in spd
        other_f_Coe = delta_y_prime / y_T  # LJ added

        f_vds=abs((y_prime[0]-y_prime[-1])/self.specimen_Data['vds'])

        g_Coe= 1 - (sum((y_prime[:-1]-y_prime[1:])**2) / sum((y_prime[:-1]-y_prime[1:]))**2 )

        q_Coe=abs(york_b)*f_Coe*g_Coe/york_sigma


        count_IZ= self.specimen_Data['steps_Arai'].count('IZ')
        count_ZI= self.specimen_Data['steps_Arai'].count('ZI')
        if count_IZ >1 and count_ZI >1:
            self.pars['magic_method_codes']="LP-PI-BT-IZZI"
        elif count_IZ <1 and count_ZI >1:
            self.pars['magic_method_codes']="LP-PI-ZI"
        elif count_IZ >1 and count_ZI <1:
            self.pars['magic_method_codes']="LP-PI-IZ"            
        else:
            self.pars['magic_method_codes']=""
            
        self.pars["specimen_b"]=york_b
        self.pars["specimen_int"]=-1*self.pars['lab_dc_field']*self.pars["specimen_b"] # possibly this is B_anc??
        self.pars["specimen_YT"]=y_T       

        self.pars["specimen_XT"] = x_T # LJ added
        self.pars['B_lab']=(float(self.specimen_Data['lab_dc_field'])) #*1e6) # LJ added.  possibly should be raw, not multiplied 1e6
        self.pars["B_anc"]= abs(york_b) * self.pars["B_lab"] #self.pars["specimen_int"] # LJ added


        self.pars["specimen_b_sigma"]=york_sigma

        self.pars["B_anc_sigma"] = york_sigma * self.pars["B_lab"]# LJ added

        self.pars["specimen_b_beta"]=beta_Coe
        self.pars["specimen_f"]=f_Coe
        self.pars["other_specimen_f"] = other_f_Coe
        self.pars["specimen_fvds"]=f_vds
        self.pars["specimen_g"]=g_Coe
        self.pars["specimen_q"]=q_Coe
        self.pars['magic_method_codes']+=":IE-TT"
        if 'x_ptrm_check' in self.specimen_Data.keys():
            if len(self.specimen_Data['x_ptrm_check'])>0:
                self.pars['magic_method_codes']+=":LP-PI-ALT-PTRM"
        if 'x_tail_check' in self.specimen_Data.keys():
            if len(self.specimen_Data['x_tail_check'])>0:
                self.pars['magic_method_codes']+=":LP-PI-BT-MD"
        print "PintPars object, self.pars after york regression: "
        print self.pars
        print "finished with York_regression()"
        print "tmin is %s, tmax is %s" %(self.tmin, self.tmax)


    def calculate_all_statistics(self):
        print "self.pars before York regression:"
        print self.pars
        print "calling calculate_all_statistics in spd.py"
        self.York_Regression()
        print "done with calculate_all_statistics"
