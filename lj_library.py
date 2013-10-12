#! /usr/bin/env python

from scipy import *
from scipy.optimize import curve_fit
#import numpy



def York_Regression(x_segment, y_segment, x_mean, y_mean, n, lab_dc_field, steps_Arai):
    x_err = x_segment - x_mean
    y_err = y_segment - y_mean
    york_b = -1* sqrt( sum(y_err**2) / sum(x_err**2) )  # averaged slope
    york_sigma= sqrt ( (2 * sum(y_err**2) - 2*york_b* sum(x_err*y_err)) / ( (n-2) * sum(x_err**2) ) )
    beta_Coe=abs(york_sigma/york_b) 
    # y_T is the intercept of the extrepolated line
    # through the center of mass (see figure 7 in Coe (1978))  
    y_T = y_mean - (york_b* x_mean)
    x_T = (-1 * y_T) / york_b  # x intercept
    # # calculate the extarplated data points for f and fvds
    x_tag=(y_segment - y_T ) / york_b # returns array of y points minus the y intercept, divided by slope
    y_tag=york_b*x_segment + y_T
    # intersect of the dashed square and the horizontal dahed line  next to delta-y-5 in figure 7, Coe (1978)
    x_prime=(x_segment+x_tag) / 2
    y_prime=(y_segment+y_tag) / 2
    delta_x_prime = abs(x_prime[-1] - x_prime[0]) #Lj add.  this is the TRM length of the best fit line
    delta_y_prime = abs(y_prime[-1] - y_prime[0]) # LJ add.  this is the NRM length of the best fit line 
    f_Coe = delta_y_prime / y_T  # LJ added 
#        g_Coe= 1 - (sum((y_prime[:-1]-y_prime[1:])**2) / sum((y_prime[:-1]-y_prime[1:]))**2 )  # old version 
    g_Coe =  1 - (sum((y_prime[:-1]-y_prime[1:])**2) / delta_y_prime ** 2)  # gap factor
    g_lim = (float(n) - 2) / (float(n) - 1)  # NOT SURE ABOUT THIS ONE 
    q_Coe=abs(york_b)*f_Coe*g_Coe/york_sigma
    w_Coe = q_Coe / sqrt(n - 2)
    count_IZ= steps_Arai.count('IZ')
    count_ZI= steps_Arai.count('ZI')
    B_lab = lab_dc_field
    B_anc = abs(york_b) * B_lab
    return {'x_err': x_err, 'y_err': y_err, 'x_tag': x_tag, 'y_tag': y_tag, 'specimen_b': york_b, 'specimen_b_sigma': york_sigma, 'specimen_b_beta': beta_Coe, 'y_int': y_T, 'x_int': x_T, 'x_prime': x_prime, 'y_prime': y_prime, 'delta_x_prime': x_prime, 'delta_y_prime': y_prime, 'specimen_f': f_Coe, 'specimen_g': g_Coe, 'specimen_g_lim': g_lim, 'specimen_q': q_Coe, 'specimen_w': w_Coe, 'count_IZ': count_IZ, 'count_ZI': count_ZI, 'B_lab': B_lab, 'B_anc': B_anc}
    

def get_vds(zdata, delta_y_prime, start, end):  # 
    """takes zdata array: [[1, 2, 3], [3, 4, 5]], delta_y_prime: 1, start value and end value.  gets vds and f_vds, etc. """
    print "calling get_vds()"
    vector_diffs = []
    for k in range(len(zdata)-1):
            # gets diff between two vectors                                                                      
        vector_diffs.append(sqrt(sum((array(zdata[k + 1 ])-array(zdata[k]))**2)))
    vector_diffs.append(sqrt(sum(array(zdata[-1])**2))) # last vector of the vds                             
    last_vector = sqrt(sum(array(zdata[-1])**2))
    vds = sum(vector_diffs)
    f_vds = abs(delta_y_prime / vds) # fvds varies, because of delta_y_prime, but vds does not.              
    vector_diffs_segment = vector_diffs[start:end]
    partial_vds = sum(vector_diffs_segment)
    max_diff = max(vector_diffs_segment)
    GAP_MAX = max_diff / partial_vds
    return {'max_diff': max_diff, 'vector_diffs': vector_diffs, 'specimen_vds': vds, 'f_vds': f_vds, 'vector_diffs_segment': vector_diffs_segment, 'partial_vds': partial_vds, 'GAP-MAX': GAP_MAX}


def get_SCAT_box(slope, slope_err, slope_beta, x_int, y_int, x_Arai_segment, y_Arai_segment, x_mean, y_mean, beta_threshold = .1):
    slope_err_threshold = abs(slope) * beta_threshold
    x, y = x_mean, y_mean
    # get lines that pass through mass center
    slope1 =  slope + (2* slope_err_threshold)
    line1_y_int = y - (slope1 * x)
    line1_x_int = -1 * (line1_y_int / slope1)
    slope2 = slope - (2 * slope_err_threshold)
    line2_y_int = y - (slope2 * x)
    line2_x_int = -1 * (line2_y_int / slope2)
        # l1_y_int and l2_x_int form the bottom line of the box       
        # l2_y_int and l1_x_int form the top line of the box              
#    print "_diagonal line1:", (0, line2_y_int), (line2_x_int, 0), (x, y)     
#    print "_diagonal line2:", (0, line1_y_int), (line1_x_int, 0), (x, y)
#    print "_bottom line:", [(0, line1_y_int), (line2_x_int, 0)]    
#    print "_top line:", [(0, line2_y_int), (line1_x_int, 0)]     
    low_bound = [(0, line1_y_int), (line2_x_int, 0)]
    high_bound = [(0, line2_y_int), (line1_x_int, 0)]
    x_max = high_bound[1][0]#              
    y_max = high_bound[0][1]
    # function for low_bound
    low_slope = (low_bound[0][1] - low_bound[1][1]) / (low_bound[0][0] - low_bound[1][0]) # y_0 - y_1 / x_0 - x_1     
    low_y_int = low_bound[0][1]
    def low_bound(x): # appears 
        y = low_slope * x + low_y_int
        return y
    # function for high_bound
    high_slope = (high_bound[0][1] - high_bound[1][1]) / (high_bound[0][0] - high_bound[1][0]) # y_0 - y_1 / x_0 - x_1     
    high_y_int = high_bound[0][1]
    def high_bound(x): 
        y = high_slope * x + high_y_int
        return y
    return low_bound, high_bound, x_max, y_max

def in_SCAT_box(x, y, low_bound, high_bound, x_max, y_max):
    """determines if a particular point falls within a box"""
    passing = True
    upper_limit = high_bound(x)
    lower_limit = low_bound(x)
    if x > x_max or y > y_max:
        print "x or y greater than x or y_max"
        passing = False
    if x < 0 or y < 0:
        print "x or y smaller than zero: all data should be positive"
        passing = False
    if y > upper_limit:
        print "y > upper limit"
        passing = False
    if y < lower_limit:
        print "y < lower_limit"
        passing = False
    return passing


def get_SCAT_points(x_Arai_segment, y_Arai_segment, tmin, tmax, ptrm_checks_temperatures, ptrm_checks_starting_temperatures, x_ptrm_check, y_ptrm_check, tail_checks_temperatures, tail_checks_starting_temperatures, x_tail_check, y_tail_check):
    """returns relevant points for a SCAT test"""
    points = []
    for i in range(len(x_Arai_segment)):
        x = x_Arai_segment[i]
        y = y_Arai_segment[i]
        points.append((x, y))

    for num, temp in enumerate(ptrm_checks_temperatures): # seems to 
        if temp >= tmin and temp <= tmax: # if temp is within selected range
            if ptrm_checks_starting_temperatures[num] >= tmin and ptrm_checks_starting_temperatures[num] <= tmax: # and also if it was not done after an out-of-range temp
                x = x_ptrm_check[num]
                y = y_ptrm_check[num]
                points.append((x, y))

    for num, temp in enumerate(tail_checks_temperatures):  # check this one   
        if temp >= tmin and temp <= tmax:
            if tail_checks_starting_temperatures[num] >= tmin and tail_checks_starting_temperatures[num] <= tmax:
                x = x_tail_check[num]
                y = y_tail_check[num]
                points.append((x, y))
           # print "points (tail checks added)", points
    return points

def get_SCAT(points, low_bound, high_bound, x_max, y_max):
    """
    runs SCAT test and returns boolean
    """
    print "hi"
    print points
    # iterate through all relevant points and see if any of them fall outside of your SCAT box
    p = True
    for point in points:
        result = in_SCAT_box(point[0], point[1], low_bound, high_bound, x_max, y_max)
        if result == False:
            print "SCAT TEST FAILED"
            p = False
    if p:
        print "SCAT TEST PASSED"
        SCAT = True
    else:
        print "SCAT TEST FAILED"
        SCAT = False
    return SCAT


def get_curve(x_Arai, y_Arai):
    
    def f(x, r, a, b):
        y = abs(sqrt(r**2-(x-a)**2)) + b
        return y
    # get best fit circle
    curve = curve_fit(f, x_Arai, y_Arai)
    r = curve[0][0] # radius of best fit circle
    a = curve[0][1] # x coordinate of circle center
    b = curve[0][2] # y coordinate of circle center
    best_fit_circle = {'a': a, 'b': b, 'radius': r}
    k = 1 /r
    # get centroid
    v = len(x_Arai)
    C_x = sum(x_Arai) / v  # x coordinate of centroid
    C_y = sum(y_Arai) / v  # y coordinate of centroid
    centroid = (C_x, C_y)
    # get direction of curvature
    if C_x < a and C_y < b:
        k = k
    if a < C_x and b < C_y:
        k = -k
    if a == C_x and b == C_y:
        k = 0
    # get SSE -- quality of best fit circle
    SSE = 0
    for i in range(len(x_Arai)):
        x = x_Arai[i]
        y = y_Arai[i]
        v = (sqrt( (x -a)**2 + (y - b)**2 ) - r )**2
#            print v                                                                                                  
        SSE += v
    return {'centroid': centroid, 'k': k, 'best_fit_circle': best_fit_circle, 'SSE': SSE }


def get_FRAC(vds, vector_diffs_segment):   
    for num in vector_diffs_segment:
        if num < 0:
            raise ValueError
    FRAC=sum(vector_diffs_segment)/ vds
    print FRAC
    return FRAC

def get_R_corr2(x_avg, y_avg, x_segment, y_segment): # 
    numerator = 0
    denominator_x = 0
    denominator_y = 0
    for num, x in enumerate(x_segment):
        r = ((x - x_avg) **2 ) * ((y_segment[num] - y_avg) **2 )
        numerator += r
#    print "numerator", numerator
    for x in x_segment:
        denominator_x += ((x - x_avg) ** 2)
#    print "den_x", denominator_x
    for y in y_segment:
        denominator_y += ((y - y_avg) ** 2)
#    print "den_y", denominator_y
    R_corr2 = numerator / (denominator_x * denominator_y)
    print R_corr2
    return R_corr2

def get_R_det2(y_segment, y_avg, y_prime):
    """
    takes in an array of y values, the mean of those values, and the array of y prime values.  returns R_det2
    """
    top = 0
    for num, y in enumerate(y_segment):
        result = (y - y_prime[num]) ** 2
        top += result
    bottom = 0
    for num, y in enumerate(y_segment):
        result = (y - y_avg) **2
        bottom += result
    print "top, bottom", top, bottom
    R_det2 = 1 - (top / bottom)
    print R_det2
    return R_det2
