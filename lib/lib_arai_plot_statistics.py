#! /usr/bin/env python

from scipy import *
import numpy


def York_Regression(x_segment, y_segment, x_mean, y_mean, n, lab_dc_field, steps_Arai):
    x_err = x_segment - x_mean
    y_err = y_segment - y_mean
    york_b = -1* sqrt( sum(y_err**2) / sum(x_err**2) )  # averaged slope

    b = sign(sum(x_err * y_err)) * std(y_segment, ddof=1)/std(x_segment, ddof=1) # ddof is degrees of freedom
    york_b = b

    york_sigma= sqrt ( (2 * sum(y_err**2) - 2*york_b* sum(x_err*y_err)) / ( (n-2) * sum(x_err**2) ) )
    beta_Coe=abs(york_sigma/york_b) 
    # y_T is the intercept of the extrepolated line
    # through the center of mass (see figure 7 in Coe (1978))  
    y_T = y_mean - (york_b* x_mean)
    x_T = (-1 * y_T) / york_b  # x intercept
    # # calculate the extrarpolated data points for f and fvds
    x_tag = (y_segment - y_T ) / york_b # returns array of y points minus the y intercept, divided by slope
    y_tag = york_b*x_segment + y_T

    # intersect of the dashed square and the horizontal dahed line  next to delta-y-5 in figure 7, Coe (1978)
    x_prime = (x_segment+x_tag) / 2  
    y_prime = (y_segment+y_tag) / 2

    delta_x_prime = abs(max(x_prime) - min(x_prime)) # TRM length of best fit line
    delta_y_prime = abs(max(y_prime) - min(y_prime)) # NRM length of best fit line

    f_Coe = delta_y_prime / abs(y_T)
    g_Coe =  1 - (sum((y_prime[:-1]-y_prime[1:])**2) / delta_y_prime ** 2)  # gap factor
    g_lim = (float(n) - 2) / (float(n) - 1) 
    q_Coe = abs(york_b)*f_Coe*g_Coe/york_sigma
    w_Coe = q_Coe / sqrt(n - 2)
    count_IZ = steps_Arai.count('IZ')
    count_ZI = steps_Arai.count('ZI')
    B_lab = lab_dc_field * 1e6
    B_anc = abs(york_b) * B_lab
    specimen_int = -1* lab_dc_field * york_b
    B_anc_sigma = york_sigma * B_lab
    return {'x_err': x_err, 'y_err': y_err, 'x_tag': x_tag, 'y_tag': y_tag, 
            'specimen_b': york_b, 'specimen_b_sigma': york_sigma, 'specimen_b_beta': beta_Coe, 
            'y_int': y_T, 'x_int': x_T, 'x_prime': x_prime, 'y_prime': y_prime, 
            'delta_x_prime': delta_x_prime, 'delta_y_prime': delta_y_prime, 'specimen_f': f_Coe, 
            'specimen_g': g_Coe, 'specimen_g_lim': g_lim, 'specimen_q': q_Coe, 'specimen_w': w_Coe, 
            'count_IZ': count_IZ, 'count_ZI': count_ZI, 'B_lab': B_lab, 'B_anc': B_anc, 
            'B_anc_sigma': B_anc_sigma, 'specimen_int': specimen_int}

def get_vds(zdata, delta_y_prime, start, end): 
    """takes zdata array: [[1, 2, 3], [3, 4, 5]], 
    delta_y_prime: 1, start value, and end value.  gets vds and f_vds, etc. """
    vector_diffs = []
    for k in range(len(zdata)-1): # gets diff between two vectors
        vector_diffs.append(sqrt( sum((array(zdata[k+1]) - array(zdata[k]))**2) ))
    last_vector = numpy.linalg.norm(zdata[-1])
    vector_diffs.append(last_vector)
    vds = sum(vector_diffs)
    f_vds = abs(delta_y_prime / vds) # fvds varies, because of delta_y_prime, but vds does not.              
    vector_diffs_segment = vector_diffs[start:end]
    partial_vds = sum(vector_diffs_segment)
    max_diff = max(vector_diffs_segment)
    GAP_MAX = max_diff / partial_vds #
    return {'max_diff': max_diff, 'vector_diffs': vector_diffs, 'specimen_vds': vds, 
            'specimen_fvds': f_vds, 'vector_diffs_segment': vector_diffs_segment, 
            'partial_vds': partial_vds, 'GAP-MAX': GAP_MAX}

def get_SCAT_box(slope,  x_mean, y_mean, beta_threshold = .1): 
    """
    takes in data and returns information about SCAT box: 
    the largest possible x_value, the largest possible y_value, 
    and functions for the two bounding lines of the box
    """
    slope_err_threshold = abs(slope) * beta_threshold
    x, y = x_mean, y_mean
    # get lines that pass through mass center, with opposite slope
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
    low_slope = (low_bound[0][1] - low_bound[1][1]) / (low_bound[0][0] - low_bound[1][0]) # 
    low_y_int = low_bound[0][1]
    def low_bound(x): 
        y = low_slope * x + low_y_int
        return y
    # function for high_bound
    high_slope = (high_bound[0][1] - high_bound[1][1]) / (high_bound[0][0] - high_bound[1][0]) # y_0-y_1/x_0-x_1     
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
        passing = False
    if x < 0 or y < 0:
        passing = False
    if y > upper_limit:
        passing = False
    if y < lower_limit:
        passing = False
    return passing

def get_SCAT_points(x_Arai_segment, y_Arai_segment, tmin, tmax, ptrm_checks_temperatures, 
                    ptrm_checks_starting_temperatures, x_ptrm_check, y_ptrm_check, 
                    tail_checks_temperatures, tail_checks_starting_temperatures, 
                    x_tail_check, y_tail_check):
    """returns relevant points for a SCAT test"""
    points = []
    for i in range(len(x_Arai_segment)): # uses only the best_fit segment, so no need for further selection
        x = x_Arai_segment[i]
        y = y_Arai_segment[i]
        points.append((x, y))

    for num, temp in enumerate(ptrm_checks_temperatures): # 
        if temp >= tmin and temp <= tmax: # if temp is within selected range
            if (ptrm_checks_starting_temperatures[num] >= tmin and 
                    ptrm_checks_starting_temperatures[num] <= tmax): # and also if it was not done after an out-of-range temp
                x = x_ptrm_check[num]
                y = y_ptrm_check[num]
                points.append((x, y))

    for num, temp in enumerate(tail_checks_temperatures):  
        if temp >= tmin and temp <= tmax:
            if (tail_checks_starting_temperatures[num] >= tmin and 
                    tail_checks_starting_temperatures[num] <= tmax):
                x = x_tail_check[num]
                y = y_tail_check[num]
                points.append((x, y))
           # print "points (tail checks added)", points
    return points

def get_SCAT(points, low_bound, high_bound, x_max, y_max):
    """
    runs SCAT test and returns boolean
    """
    # iterate through all relevant points and see if any of them fall outside of your SCAT box
    SCAT = True
    for point in points:
        result = in_SCAT_box(point[0], point[1], low_bound, high_bound, x_max, y_max)
        if result == False:
           # print "SCAT TEST FAILED"
            SCAT = False
    return SCAT

def get_FRAC(vds, vector_diffs_segment):   
    for num in vector_diffs_segment:
        if num < 0:
            raise ValueError('vector diffs should not be negative')
    if vds == 0:
        raise ValueError('attempting to divide by zero. vds should be a positive number')
    FRAC = sum(vector_diffs_segment) / vds 
    return FRAC

def get_R_corr2(x_avg, y_avg, x_segment, y_segment): # 
    xd = x_segment - x_avg # detrend x_segment
    yd = y_segment - y_avg # detrend y_segment
    if sum(xd**2) * sum(yd**2) == 0:
        raise ValueError('attempting to divide by zero')
    rcorr = sum((xd * yd))**2 / (sum(xd**2) * sum(yd**2))
    print 'xd', xd
    print 'yd', yd
    return rcorr

def get_R_det2(y_segment, y_avg, y_prime):
    """
    takes in an array of y values, the mean of those values, and the array of y prime values.  
    returns R_det2
    """
    numerator = sum((numpy.array(y_segment) - numpy.array(y_prime))**2)
    denominator = sum((numpy.array(y_segment) - y_avg)**2)
    R_det2 = 1 - (numerator / denominator)
    return R_det2

def get_b_wiggle(x, y, y_int):
    """returns instantaneous slope from the ratio of NRM lost to TRM gained at the ith step"""
    if x == 0:
        b_wiggle = 0
    else:
        b_wiggle = (y_int - y) / x
    return b_wiggle

def get_Z(x_segment, y_segment, x_int, y_int, slope):
    """Arai plot zigzag parameter"""
    Z = 0
    first_time = True
    for num, x in enumerate(x_segment):
        b_wiggle = get_b_wiggle(x, y_segment[num], y_int)
        z = (x * abs(b_wiggle - abs(slope)) ) / abs(x_int)
        Z += z
        first_time = False
    return Z

def get_Zstar(x_segment, y_segment, x_int, y_int, slope, n):
    """Arai plot zigzag parameter (alternate)"""
    total = 0
    first_time = True
    for num, x in enumerate(x_segment):
        b_wiggle = get_b_wiggle(x, y_segment[num], y_int)
        result = 100 * ( (x * abs(b_wiggle - abs(slope)) ) / abs(y_int) )
        total += result
        first_time = False
    Zstar = (1. / (n - 1.)) * total
    return Zstar


# IZZI_MD (mainly)

def get_normed_points(point_array, norm): # good to go
#    takes a set of points and norms them
    norm = float(norm)
    floated_array = []
    for p in point_array: # need to make sure each point is a float
        floated_array.append(float(p))
    points = numpy.array(floated_array) / norm
    return points

def get_xy_array(x_segment, y_segment): 
 #takes lists of x and y coordiantes and combines them, returning: [(x, y), (x, y)]
    xy_array = []
    for num, x in enumerate(x_segment):
        xy_array.append((x, y_segment[num]))
    return xy_array

