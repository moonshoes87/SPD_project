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
        print self.s
        self.datablock = self.specimen_Data['datablock']

        self.x_Arai=self.specimen_Data['x_Arai']
        self.y_Arai=self.specimen_Data['y_Arai']
        self.t_Arai=self.specimen_Data['t_Arai']

        self.zdata = self.specimen_Data['zdata'] # LJ add

        self.x_tail_check=self.specimen_Data['x_tail_check']
        self.y_tail_check=self.specimen_Data['y_tail_check']

        self.zijdblock=self.specimen_Data['zijdblock']        
        self.z_temperatures=self.specimen_Data['z_temp']

        # tmax, tmin, start, and end -- inclusive, or exclusive???
        self.start=self.t_Arai.index(tmin)
        self.end=self.t_Arai.index(tmax)

        self.length = abs(self.end - self.start) + 1  # plus one because of indexing -- starts at 0

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
        self.stuff = ["s", "datablock", "x_Arai", "y_Arai", "t_Arai", "x_tail_check", "y_tail_check", "zijdblock", "z_temperatures", "start", "end", "pars", "specimen_Data", "tmin", "tmax", "tmin_K", "tmax_K"] # needs to be updated
 
        #LJ ADDING stats:
        self.n_max = len(self.t_Arai)  # gets total number of temperatures taken.  (p. 4, top)
        self.tmin = tmin # self-explanatory
        self.tmax = tmax
        self.tmin_K = tmin - 273.15
        self.tmax_K = tmax - 273.15

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
        y_T = y_Arai_mean - (york_b* x_Arai_mean)

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

  #      f_Coe=abs((y_prime[0]-y_prime[-1])/y_T)  # same as 'f' in spd
        f_Coe = delta_y_prime / y_T  # LJ added

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
        self.pars["specimen_int"]=-1*self.pars['lab_dc_field']*self.pars["specimen_b"] # same as B_anc
        self.pars["specimen_YT"]=y_T       

        self.pars["specimen_XT"] = x_T # LJ added
        self.pars['B_lab']=(float(self.specimen_Data['lab_dc_field'])) #*1e6) # LJ added.  possibly should be raw, not multiplied 1e6
        self.pars["B_anc"]= abs(york_b) * self.pars["B_lab"] #self.pars["specimen_int"] # LJ added


        self.pars["specimen_b_sigma"]=york_sigma

        self.pars["B_anc_sigma"] = york_sigma * self.pars["B_lab"]# LJ added

        self.pars["specimen_b_beta"]=beta_Coe
        self.pars["specimen_f"]=f_Coe

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


    def get_vds(self):  # appears now to be working.  fetches vector difference sum.  unless it needs to use start and end instead of using all the points
        """gets vds and f_vds"""
        print "calling new_get_vds()"
        zdata = self.specimen_Data['zdata']
        vector_diffs = []
        for k in range(len(zdata)-1): # (self.length -1): # but should it be this??
          #  vector_diff = (sqrt(sum((array(zdata[k + 1])-array(zdata[k]))**2)))  
           # print "vector diff: ", vector_diff
            vector_diffs.append(sqrt(sum((array(zdata[k + 1 ])-array(zdata[k]))**2)))
        vector_diffs.append(sqrt(sum(array(zdata[-1])**2))) # last vector of the vds
        last_vector = sqrt(sum(array(zdata[-1])**2)) 
        vds = sum(vector_diffs)
        print "correct vds: ", self.specimen_Data['vds']
        print "vds", vds
        delta_y_prime = self.pars['delta_y_prime']   # this is definitely correct
        f_vds = abs(delta_y_prime / vds) # fvds varies, because of delta_y_prime, but vds does not.  

# should be either vector_diffs is always all inclusive, and then you pick your range: vector_diffs[1:10], OR you calculate the VDS using start and end and then never need a vector_diffs_segment.  email Ron about this.

        self.pars['vector_diffs'] = vector_diffs
        self.pars['specimen_vds'] = vds
        self.pars["specimen_fvds"]=f_vds 

    def get_FRAC(self):   # seems to work
        vds = self.pars['specimen_vds']
        vector_diffs = self.pars['vector_diffs']
        vector_diffs_segment = vector_diffs[self.start:self.end]
        print "vector diffs_segment:", vector_diffs_segment
        FRAC=sum(vector_diffs_segment)/ vds
        print "FRAC: ", FRAC
        self.pars['FRAC'] = FRAC

# stolen from thellier_gui:
#        vector_diffs=self.Data[s]['vector_diffs']
#        vector_diffs_segment=vector_diffs[zstart:zend]
#        FRAC=sum(vector_diffs_segment)/self.Data[s]['vds']
#        max_FRAC_gap=max(vector_diffs_segment/sum(vector_diffs_segment))



    def calculate_all_statistics(self):
        print "self.pars before York regression:"
        print self.pars
        print "calling calculate_all_statistics in spd.py"
        self.York_Regression()
        self.get_vds()
        self.get_FRAC()
        print "done with calculate_all_statistics"


# C temps: [273, 373.0, 423.0, 473.0, 498.0, 523.0, 548.0, 573.0, 598.0, 623.0, 648.0, 673.0, 698.0, 723.0, 748.0, 773.0, 798.0, 823.0]
if True:
    import new_lj_thellier_gui_spd as tgs
    gui = tgs.Arai_GUI()
    thing = PintPars(gui.Data, '0238x5721063', 273., 823.)
    thing.calculate_all_statistics()
    thing.get_FRAC()
    thing1 = PintPars(gui.Data, '0238x5721063', 598., 698.)
    thing1.calculate_all_statistics()
    thing1.get_FRAC()
    print "thing f_vds: ", thing.pars['specimen_fvds']
    print thing.specimen_Data['vds']
    print "thing1 f_vds: ", thing1.pars['specimen_fvds']
    print thing1.specimen_Data['vds']
    print "thing FRAC: ", thing.pars['FRAC']
    print "thing1 FRAC: ", thing1.pars['FRAC']
    print "thing1 temps: ", thing1.tmin_K, thing1.tmax_K



if False:
    # old code that uses zijdblock instead of zdata to get vds/fvds
    def old_get_vds(self): # stolen from tgs
        """vector difference sum"""
        zijdblock = self.zijdblock
        print "zijdblock", zijdblock
        z_temperatures=[row[0] for row in zijdblock]
        print "z_temperatures", z_temperatures
        zdata=[]
        vector_diffs=[]
        NRM=zijdblock[0][3]
        print "NRM", NRM
# each zijdblock list item:  [treatment (temperature),declination,inclination,intensity,ZI,measurement_flag,magic_instrument_codes]
        for k in range(len(zijdblock)):
            DIR=[zijdblock[k][1],zijdblock[k][2],zijdblock[k][3]/NRM]  # translation dec, inc, magnetic measurement into x, y, z
           # print "treatment", zijdblock[k][0]
  #          print "ZI", zijdblock[k][4]
   #         print "DIR (vectors). dec, inc, intensity ", DIR
            cart=self.dir2cart(DIR)
#            print "cartesian coordinates, " , cart
            zdata.append(array([cart[0],cart[1],cart[2]]))
            if k>0:
                vector_diff = (sqrt(sum((array(zdata[-2])-array(zdata[-1]))**2)))
                print "vector_diff: ", vector_diff
                vector_diffs.append(sqrt(sum((array(zdata[-2])-array(zdata[-1]))**2)))
        vector_diffs.append(sqrt(sum(array(zdata[-1])**2))) # last vector of th
        last_vector = sqrt(sum(array(zdata[-1])**2)) # last vector of th
        print "last_vector", last_vector
        vds=sum(vector_diffs)  # vds calculation                    
        print "vector_diffs", vector_diffs
        print "correct vds", self.specimen_Data['vds']
        print "vds ", vds
#        print "zdata[0]", zdata[0]
        zdata=array(zdata)
 #       print "zdata", zdata

        delta_y_prime = self.pars['delta_y_prime']   # this is definitely correct
        f_vds=abs( delta_y_prime / vds)
#        vds = self.specimen_Data['vds'] # LJ added
#        self.pars['vds'] = vds # lj added
        self.pars['specimen_vds_new'] = vds
#        f_vds=abs((y_prime[0]-y_prime[-1])/self.specimen_Data['vds']) # old code, not encapsulated in spd.py
 #       self.pars["specimen_fvds"]=f_vds
        self.pars["specimen_fvds_new"]=f_vds 
        return vector_diffs, vds



    def dir2cart(self, d):
#        print "calling dir2cart(), not in anything"                                                             
       # converts list or array of vector directions, in degrees, to array of cartesian coordinates, in x,y,z    
        ints=ones(len(d)).transpose() # get an array of ones to plug into dec,inc pairs                
        d=array(d)
        rad=pi/180.
        if len(d.shape)>1: # array of vectors
            decs,incs=d[:,0]*rad,d[:,1]*rad
            if d.shape[1]==3: ints=d[:,2] # take the given lengths
        else: # single vector
            decs,incs=array(d[0])*rad,array(d[1])*rad
            if len(d)==3:
                ints=array(d[2])
            else:
                ints=array([1.])
        cart= array([ints*cos(decs)*cos(incs),ints*sin(decs)*cos(incs),ints*sin(incs)]).transpose()
        return cart
    
