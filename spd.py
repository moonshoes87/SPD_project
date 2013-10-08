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

        self.zijdblock=self.specimen_Data['zijdblock']        
        self.z_temperatures=self.specimen_Data['z_temp']

        # tmax, tmin, start, and end -- inclusive, or exclusive???
        self.start=self.t_Arai.index(tmin)
        self.end=self.t_Arai.index(tmax)

        self.length = abs(self.end - self.start) + 1  # plus one because of indexing -- starts at 0

        self.pars={}

        self.pars['lab_dc_field']=self.specimen_Data['pars']['lab_dc_field']
  #      self.pars['magic_method_codes']=Data[self.s]['pars']['magic_method_codes']
        # for some reason missing any magic_method_codes.  possibly these would have been incorporated into the data from rmag_anisotropy or something
        # magic_method codes are locked up in datablock, not actually extracted.  not sure if this happens somewhere else in thellier_gui or not
        # also, fix the weirdness of having to set the precise number for tmin and tmax
        self.pars['specimen_int_n']=self.end-self.start+1
        self.stuff = ["s", "datablock", "x_Arai", "y_Arai", "t_Arai", "x_tail_check", "y_tail_check", "tail_checks_temperatures", "tail_checks_starting_temperatures", "x_ptrm_check", "y_ptrm_check", "ptrm_checks_temperatures", "ptrm_checks_starting_temperatures", "zijdblock", "z_temperatures", "start", "end", "pars", "specimen_Data", "tmin", "tmax", "tmin_K", "tmax_K"] # needs to be updated
 
        #LJ ADDING stats:
        self.n = float(self.end-self.start+1)  # n is the number of points involved
        self.n_max = len(self.t_Arai)  # gets total number of temperatures taken.  (p. 4, top)
        self.tmin = tmin # self-explanatory
        self.tmax = tmax
        self.tmin_K = tmin - 273.15
        self.tmax_K = tmax - 273.15
        self.x_Arai_segment = self.x_Arai[self.start:self.end+1]  # returns array of relevant x points
        self.y_Arai_segment = self.y_Arai[self.start:self.end+1]


    # eventually add a function to change tmax and/or tmin.  must also change start, end, 

    def York_Regression(self):
        #-------------------------------------------------
        # York regresssion (York, 1967) following Coe et al.(1978)
        # calculate statsitics: n,b,b_sigma_,f,fvds,specimen_int,g,q
        # modified from pmag.py
        #-------------------------------------------------               

        x_Arai_segment= self.x_Arai_segment #self.x_Arai[self.start:self.end+1]  # returns array of relevant x points
        y_Arai_segment= self.y_Arai_segment #self.y_Arai[self.start:self.end+1]

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

#        g_Coe= 1 - (sum((y_prime[:-1]-y_prime[1:])**2) / sum((y_prime[:-1]-y_prime[1:]))**2 )  # old version
        g_Coe =  1 - (sum((y_prime[:-1]-y_prime[1:])**2) / delta_y_prime ** 2)  # gap factor
        g_lim = (float(n) - 2) / (float(n) - 1)  # NOT SURE ABOUT THIS ONE

        q_Coe=abs(york_b)*f_Coe*g_Coe/york_sigma

        w_Coe = q_Coe / sqrt(self.n - 2)

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
        self.pars["specimen_g_lim"] = g_lim
        self.pars["specimen_q"]=q_Coe
        self.pars["specimen_w"]= w_Coe # LJ add. untested, but it should work

        

        self.pars['magic_method_codes']+=":IE-TT"
        if 'x_ptrm_check' in self.specimen_Data.keys():
            if len(self.specimen_Data['x_ptrm_check'])>0:
                self.pars['magic_method_codes']+=":LP-PI-ALT-PTRM"
        if 'x_tail_check' in self.specimen_Data.keys():
            if len(self.specimen_Data['x_tail_check'])>0:
                self.pars['magic_method_codes']+=":LP-PI-BT-MD"
#        print "PintPars object, self.pars after york regression: "
#        print self.pars
#        print "finished with York_regression()"
#        print "tmin is %s, tmax is %s" %(self.tmin, self.tmax)


    def get_vds(self):  # appears now to be working.  fetches vector difference sum.  vds and fvds are correct.  
        """gets vds and f_vds"""
        print "calling get_vds()"
        zdata = self.specimen_Data['zdata']
        vector_diffs = []
        for k in range(len(zdata)-1): 
            # gets diff between two vectors
            vector_diffs.append(sqrt(sum((array(zdata[k + 1 ])-array(zdata[k]))**2)))
        vector_diffs.append(sqrt(sum(array(zdata[-1])**2))) # last vector of the vds
        last_vector = sqrt(sum(array(zdata[-1])**2)) 
        vds = sum(vector_diffs)
        delta_y_prime = self.pars['delta_y_prime']  
        f_vds = abs(delta_y_prime / vds) # fvds varies, because of delta_y_prime, but vds does not.  
        vector_diffs_segment = vector_diffs[self.start:self.end]
        partial_vds = sum(vector_diffs_segment)
        max_diff = max(vector_diffs_segment)
        GAP_MAX = max_diff / partial_vds
        self.pars['max_diff'] = max_diff
        self.pars['vector_diffs'] = vector_diffs
        self.pars['specimen_vds'] = vds
        self.pars["specimen_fvds"]=f_vds 
        self.pars['vector_diffs_segment'] = vector_diffs_segment
        self.pars['partial_vds'] = partial_vds
        self.pars['GAP-MAX'] = GAP_MAX


    def get_FRAC(self):   # works
        vds = self.pars['specimen_vds']
#        vector_diffs = self.pars['vector_diffs']
#        vector_diffs_segment = vector_diffs[self.start:self.end]
        vector_diffs_segment = self.pars['vector_diffs_segment']  # could use this instead
        FRAC=sum(vector_diffs_segment)/ vds
        self.pars['FRAC'] = FRAC


    def get_curve(self):  # need to check this shit with Alex.  also, see about somewebsite.com/code
        # a, b == x, y coordinates of the circle center
        # x is the x coordinate of some point
        # so x coordinate of the circle center * hyperbolic tangent of (y coordinate of center * x)
#        def tan_h(x, a, b):
#            return a*tanh(b*x)

#        def scipy.optimize.curve_fit(f, xdata, ydata)
        def f(x, r, a, b): # circle function
            y = abs(sqrt(r**2-(x-a)**2)) + b
            return y
        import scipy
        import numpy
        curve = scipy.optimize.curve_fit(f, self.x_Arai, self.y_Arai)
        r = curve[0][0] # radius of circle
        a = curve[0][1] # x coordinate of circle center
        b = curve[0][2] # y coordinate of circle center
        k = 1/r
        # get data centroid
#        centroid = []
#        for n in range(len(self.x_Arai)): # possibly this needs to be a smaller segment, i.e. using self.start:self.end
#            point = [self.x_Arai[n], self.y_Arai[n]]
#           # print "x", x
#            centroid.append(point)
#        centroid = numpy.array(centroid)   
#        print "centroid of data points:", centroid
        v = len(self.x_Arai)  # would have to change this to reflect the proper length, also
#        centroid_x_sum = centroid[:,0].sum()  # sums x values
#        centroid_y_sum = centroid[:,1].sum()  # sums y values
#        centroid = numpy.array([centroid_x_sum, centroid_y_sum]) / v  # divides by the number of data points to find the "average" poitn
#        C_x = centroid[0] # x coordinate of centroid
#        C_y = centroid[1] # y coordinate of centroid
        # possibly simpler:  # but this too might need self.start:self.end
        C_x = sum(self.x_Arai) / v # x coordinate of centroid
        C_y = sum(self.y_Arai) / v # y coordinate of centroid
        # done getting centroid
        # get "direction" of the curvature
        if C_x < a and C_y < b:
            k = k
        if a < C_x and b < C_y:
            k = -k
        if a == C_x and b == C_y:
            k = 0
        SSE = 0 # quality of best_fit circle
        for i in range(len(self.x_Arai)):
            x = self.x_Arai[i]
            y = self.y_Arai[i]
            v = (sqrt( (x -a)**2 + (y - b)**2 ) - r )**2
#            print v
            SSE += v
#        print SSE
        self.pars['centroid'] = (C_x, C_y)
        self.pars['specimen_k'] = k
        self.pars['best_fit_circle'] = { "a": a, "b" : b, "radius": r }
        self.pars['SSE'] = SSE
        return k, a, b

    def get_SCAT(self):
        slope, slope_err, beta = self.pars['specimen_b'], self.pars['specimen_b_sigma'], self.pars['specimen_b_beta']
# need beta_threshold.  default in thellier_gui is .1
        x_int, y_int = self.pars['specimen_XT'], self.pars['specimen_YT']
        beta_threshold = 0.1 
        slope_err_threshold = abs(slope) * beta_threshold
        x_mean=mean(array(self.x_Arai_segment))
        y_mean=mean(array(self.y_Arai_segment))
        mass_center = (x_mean, y_mean)
        print "center of mass:", mass_center # 
        # mass center CAN be different from centroid previously calculated, because SCAT allows for a smaller subset and best fit circle or whatever does not
        # possibly the above is incorrect.  possibly you should readjust centroid in get_curve to allow for a subset of the data.  not sure.

        # getting lines passing through mass_center
        x = mass_center[0]
        y = mass_center[1]
#        print "beta_threshold", beta_threshold
#        print "slope_err_threshold", slope_err_threshold
        slope_1 = slope + (2 * slope_err_threshold)
        l1_y_int = y - (slope_1 * x)
        l1_x_int = -1 * (l1_y_int / slope_1)
        slope_2 = slope - 2 * slope_err_threshold
        l2_y_int = y - (slope_2 * x)
        l2_x_int = -1 * (l2_y_int / slope_2)
        # l1_y_int and l2_x_int form the bottom line of the box
        # l2_y_int and l1_x_int form the top line of the box
#        print "diagonal line1:", (0, l2_y_int), (l2_x_int, 0), (x, y)
#        print "diagonal line2:", (0, l1_y_int), (l1_x_int, 0), (x, y)
#        print "center of mass: ", mass_center
#        print "bottom line:", [(0, l1_y_int), (l2_x_int, 0)]
#        print "top line:", [(0, l2_y_int), (l1_x_int, 0)]
        low_line = [(0, l1_y_int), (l2_x_int, 0)]
        high_line = [(0, l2_y_int), (l1_x_int, 0)]
        x_max = high_line[1][0]# 
        y_max = high_line[0][1]
        print "Max x, y: ", x_max, y_max
        
#        function for bottom line:
        low_slope = (low_line[0][1] - low_line[1][1]) / (low_line[0][0] - low_line[1][0]) # y_0 - y_1 / x_0 - x_1
        low_y_int = low_line[0][1]
#        y = mx + b
        def low_line(x): # appears correct
            y = low_slope * x + low_y_int
            return y

        # function for top line
        high_slope = (high_line[0][1] - high_line[1][1]) / (high_line[0][0] - high_line[1][0]) # y_0 - y_1 / x_0 - x_1
        high_y_int = high_line[0][1]
        def high_line(x): # appears correct
            y = high_slope * x + high_y_int
            return y
#        print "High:"
#        print "x = 1", high_line(1.)
#        print "x =2", high_line(2.)
#        print "x=3", high_line(3.)
#        print "low"
#        print "x=.5", low_line(.5)
#        print "x=1.5", low_line(1.5)
#        print "x=3", low_line(3.)
#        print "-"

        def in_SCAT(x, y):
            print "x, y", x, y
            passing = True
            if x > x_max or y > y_max: 
                print "x or y greater than x or y_max"
                passing = False
            if x < 0 or y < 0: # all data must be in upper right quadrant of graph
                print "x or y smaller than 0"
                passing = False
            upper_limit = high_line(x)
            lower_limit = low_line(x)
            print "upper limit, lower limit for y: ", upper_limit, lower_limit
            if y > upper_limit:
                print "y > upper limit"
                passing = False
            if y < lower_limit:
                print "y < lower limit"
                passing = False
           # print "SCAT pass is: ", passing
            return passing # boolean
    

#        line1 -- passes through mass center with slope of: slope + 2(beta_threshold)
        # line 2 -- passes through mass center with slope of: slope - 2(beta_threshold)
        # forms SCAT box with x and y intercepts
        
        # first: figure out which points to include (slightly complex regarding pTRM checks)
        # figure out how to determine if they are "within the box"
        # if all are in box, then SCAT is "TRUE"
        points = []
        # add these together:
        
#        for i in range(len(self.x_Arai[self.start:self.end+1])):
        for i in range(self.start, self.end + 1):
            x = self.x_Arai[i]
            y = self.y_Arai[i]
            points.append((x, y))
        
        # for loop that appends proper, relevant ptrm checks

        print "points (x and y Arai added)", points

        for num, temp in enumerate(self.ptrm_checks_temperatures): # seems to work
            if temp >= self.tmin and temp <= self.tmax: # if temp is within selected range
                if self.ptrm_checks_starting_temperatures[num] >= self.tmin and self.ptrm_checks_starting_temperatures[num] <= self.tmax: # and also if it was not done after an out-of-range temperature
                    x = self.x_ptrm_check[num]
                    y = self.y_ptrm_check[num]
                    points.append((x, y))

        print "points (ptrm checks added)", points
        # for loop that appends proper, relevant tail checks
        for num, temp in enumerate(self.tail_checks_temperatures):  # check this one
            if temp >= self.tmin and temp <= self.tmax:
                if self.tail_checks_starting_temperatures[num] >= self.tmin and self.tail_checks_starting_temperatures[num] <= self.tmax:
                    x = self.x_tail_check[num]
                    y = self.y_tail_check[num]
                    points.append((x, y))
        print "points (tail checks added)", points

        # iterate through points and see if any of them are outside of your SCAT box
        p = True
        for point in points:
            result = in_SCAT(point[0], point[1])
            if result == False:
                print "SCAT TEST FAILED"
                p = False
        if p:
            print "SCAT TEST PASSED"
        else:
            print "SCAT TEST FAILED"
        print "---------"
        
        
                            

# pTRM checks ("triangles")
#         Data[s]['x_ptrm_check']=[] # a list of x coordinates of pTRM checks
#         Data[s]['y_ptrm_check']=[] # a list of y coordinates of pTRM checks      
#         Data[s]['ptrm_checks_temperatures']=[] # a list of pTRM checks temperature 
#         Data[s]['x_ptrm_check_starting_point'] # a list of x coordinates of the point ehere the pTRM checks started from
#         Data[s]['y_ptrm_check_starting_point'] # a list of y coordinates of the point ehere the pTRM checks started from            #         Data[s]['ptrm_checks_starting_temperatures']=[] # a list of temperatures from which the pTRM checks started from 


        # for loop that appends ptrm tail checks
        

        # make a function that puts together x_Arai and y_Arai points, since you use that multiple times

            
        
        
    def calculate_all_statistics(self):
        print "calling calculate_all_statistics in spd.py"
        self.York_Regression()
        self.get_vds()
        self.get_FRAC()
        self.get_curve()
        print "done with calculate_all_statistics"

# K temps: [0.0, 100.0, 150.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0]
# C temps: [273, 373.0, 423.0, 473.0, 498.0, 523.0, 548.0, 573.0, 598.0, 623.0, 648.0, 673.0, 698.0, 723.0, 748.0, 773.0, 798.0, 823.0]
if True:
    import new_lj_thellier_gui_spd as tgs
    gui = tgs.Arai_GUI()
    specimens = gui.Data.keys()
    thing = PintPars(gui.Data, '0238x6011044', 473., 648.)  # this one only gets incorrect results
    thing.calculate_all_statistics()
#    thing1 = PintPars(gui.Data, specimens[3], 523., 773.)
#    thing1.calculate_all_statistics()
#    thing2 = PintPars(gui.Data, specimens[4], 273., 798.)
#    thing2.calculate_all_statistics()
#    thing3 = PintPars(gui.Data, specimens[5], 598, 698)
#    thing3.calculate_all_statistics()
    print "---"
    print thing.s, thing.tmin_K, thing.tmax_K  # incorrect
    thing.get_SCAT()
#    print thing1.s, thing1.tmin_K, thing1.tmax_K
#    thing1.get_SCAT()
#    print thing2.s, thing2.tmin_K, thing2.tmax_K
#    thing2.get_SCAT()
#    print thing3.s, thing3.tmin_K, thing3.tmax_K
#    thing3.get_SCAT()
    
    


##    print "thing f_vds: ", thing.pars['specimen_fvds']
 #   print thing.specimen_Data['vds']
#    print "thing1 f_vds: ", thing1.pars['specimen_fvds']
#    print thing1.specimen_Data['vds']
 #   print "thing1 temps: ", thing1.tmin_K, thing1.tmax_K
  #  print "thing2 f_vds: ", thing2.pars['specimen_fvds']
#    print "thing2.s: ", thing2.s
#    print "thing2 tmin and tmax: ", thing2.tmin_K, thing2.tmax_K
#    print "thing2 GAP-MAX: ", thing2.pars['alt_GAP-MAX'] 
 #   print "thing3 f_vds: ", thing3.pars['specimen_fvds']
#    print "thing3.s: ", thing3.s
#    print "thing3 tmin and tmax: ", thing3.tmin_K, thing3.tmax_K
#    print "thing3 GAP-MAX: ", thing3.pars['alt_GAP-MAX'] 
 




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
    
