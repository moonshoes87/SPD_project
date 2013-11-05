#! /usr/bin/env python

from scipy import *
from scipy.optimize import curve_fit
import numpy



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
    g_lim = (float(n) - 2) / (float(n) - 1) 
    q_Coe=abs(york_b)*f_Coe*g_Coe/york_sigma
    w_Coe = q_Coe / sqrt(n - 2)
    count_IZ= steps_Arai.count('IZ')
    count_ZI= steps_Arai.count('ZI')
    B_lab = lab_dc_field
    B_anc = abs(york_b) * B_lab
    specimen_int = -1* lab_dc_field * york_b
    B_anc_sigma = york_sigma * B_lab
    return {'x_err': x_err, 'y_err': y_err, 'x_tag': x_tag, 'y_tag': y_tag, 'specimen_b': york_b, 'specimen_b_sigma': york_sigma, 'specimen_b_beta': beta_Coe, 'y_int': y_T, 'x_int': x_T, 'x_prime': x_prime, 'y_prime': y_prime, 'delta_x_prime': delta_x_prime, 'delta_y_prime': delta_y_prime, 'specimen_f': f_Coe, 'specimen_g': g_Coe, 'specimen_g_lim': g_lim, 'specimen_q': q_Coe, 'specimen_w': w_Coe, 'count_IZ': count_IZ, 'count_ZI': count_ZI, 'B_lab': B_lab, 'B_anc': B_anc, 'B_anc_sigma': B_anc_sigma, 'specimen_int': specimen_int}
    

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
    return {'max_diff': max_diff, 'vector_diffs': vector_diffs, 'specimen_vds': vds, 'specimen_fvds': f_vds, 'vector_diffs_segment': vector_diffs_segment, 'partial_vds': partial_vds, 'GAP-MAX': GAP_MAX}


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

    for num, temp in enumerate(tail_checks_temperatures):  
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
            raise ValueError('vector diffs should not be negative')
    if vds == 0:
        raise ValueError('attempting to divide by zero. vds should be a positive number')
    FRAC=sum(vector_diffs_segment)/ vds
   # print FRAC
    return FRAC

def get_R_corr2(x_avg, y_avg, x_segment, y_segment): # 
    numerator = 0
    denominator_x = 0
    denominator_y = 0
    for num, x in enumerate(x_segment):
        r = ((x - x_avg)  * (y_segment[num] - y_avg) )
        numerator += r
#    print "numerator", numerator
    numerator = numerator**2
    for x in x_segment:
        denominator_x += ((x - x_avg) ** 2)
#    print "den_x", denominator_x
    for y in y_segment:
        denominator_y += ((y - y_avg) ** 2)
#    print "den_y", denominator_y
    denominator = denominator_x * denominator_y
    if denominator == 0: raise ValueError("get_R_corr2 attempted dividing by zero")
    R_corr2 = numerator / denominator
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
    R_det2 = 1 - (top / bottom)
    return R_det2

def get_b_wiggle(x, y, y_int):
    """returns instantaneous slope from the ratio of NRM lost to TRM gained at the ith step"""
#    print "(y_int - y) / x", 
#    print "(", y_int,  "-", y, ") /",  x
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
        z = (x * abs(b_wiggle - abs(slope)) ) / x_int
        Z += z
        first_time = False
    return Z

def get_Zstar(x_segment, y_segment, x_int, y_int, slope, n):
    """Arai plot zigzag parameter (alternate)"""
#    for num, x in enumerate(x_segment):
#        b_wiggle = get_b_wiggle(x, y_segment[num], y_int)
    print "doing Zstar"
    total = 0
    first_time = True
    for num, x in enumerate(x_segment):
        b_wiggle = get_b_wiggle(x, y_segment[num], y_int)
        result = 100 * ( (x * abs(b_wiggle - abs(slope)) ) / y_int )
        total += result
        first_time = False
    Zstar = (1. / (n - 1.)) * total
    print "end Zstar"
    return Zstar

# IZZI_MD (mainly)

def get_normed_points(point_array, norm): # good to go
    """takes a set of points and norms them"""
    norm = float(norm)
    floated_array = []
    for p in point_array: # need to make sure each point is a float
        floated_array.append(float(p))
    points = numpy.array(floated_array) / norm
    return points

def get_xy_array(x_segment, y_segment): 
    """takes lists of x and y coordiantes and combines them, returning: [(x, y), (x, y)]"""
    xy_array = []
    for num, x in enumerate(x_segment):
        xy_array.append((x, y_segment[num]))
    return xy_array



xy_segment = [(1, 2), (3, 4), (5, 6), (7, 8), (9,10), (11, 12), (13, 14)]
steps = ['IZ', 'ZI', 'ZI', 'IZ', 'IZ', 'ZI', 'IZ']

def get_triangles(xy_segment = xy_segment, Arai_steps = steps): # works!
    segment = xy_segment[1: ]
    no_repeat_segment = []
    no_repeat_steps = []
    last_step = ""
    for i in segment:
        index = xy_segment.index(i) # must be xy_segment to correspond to the Arai_steps, not the truncated segment
        step = Arai_steps[index]
        if step == last_step:
            no_repeat_segment[-1] = i
            no_repeat_steps[-1] = step
        else:
            no_repeat_segment.append(i)
            no_repeat_steps.append(step)
        last_step = step
    triangles = []
    midpoints = []
    for num, i in enumerate(no_repeat_segment[:-2]):
        triangles.append((i, no_repeat_segment[num + 1], no_repeat_segment[num + 2]))
        midpoints.append(no_repeat_steps[num+1])
    print "triangles:", triangles
    return {'triangles': triangles, 'midpoints': midpoints}

    
def get_triangle_sides(x_segment, y_segment):
    """finds the length of the sides of a triangle from three sets of x, y coordinates"""
    L1 = sqrt((x_segment[0] - x_segment[1])**2 + (y_segment[0] - y_segment[1])**2)
    L2 = sqrt((x_segment[1] - x_segment[2])**2 + (y_segment[1] - y_segment[2])**2)
    L3 = sqrt((x_segment[2] - x_segment[0])**2 + (y_segment[2] - y_segment[0])**2)
    return {'L1': L1, 'L2': L2, 'L3': L3}


def get_triangle(line1, line2, line3):
    """takes length of a triangle's lines and returns angle1, the triangle's height, and its area"""
#    print "top of phi: ", (line2**2 + line3**2 - line1**2)
#    print "bottom of phi: ", 2 * line2 * line3
#    print "arccos of all that"
    phi = arccos((line2**2 + line3**2 - line1**2) / (2 * line2 * line3))
    height = line3 * sin(phi)
    area = (line2 * line3 * sin(phi)) / 2
    return { 'triangle_phi': phi, 'triangle_H': height, 'triangle_A': area }
       

def get_sign(triangle, midpoint):
    first_line = [(triangle[0]), (triangle[2])]
    first_slope = (first_line[1][1] - first_line[0][1]) / (first_line[1][0] - first_line[0][0])
    first_y_int = first_line[0][1] - (first_slope * first_line[0][0])
    second_line = [(triangle[1])]
    second_y_int = second_line[0][1] - (first_slope * second_line[0][0])
    sign = 0.
    if midpoint == 'IZ':
        if first_y_int > second_y_int:
            sign = 1.
        elif first_y_int < second_y_int:
            sign = -1.
        else:
            sign = 0
    if midpoint == 'ZI':
        if first_y_int > second_y_int:
            sign = -1.
        elif first_y_int < second_y_int:
            sign = 1.
        else:
            sign = 0
#    return sign
    return {'slope': first_slope, 'first_y_int': first_y_int, 'second_y_int': second_y_int, 'sign': sign }


x = numpy.array([1., 1.5, 2., 3., 6., 8., 8.5, 8.5])
y = numpy.array([8., 4., 3.9, 3.5, 3., 2.5, 2., 1.5])
xy = []
for num, x in enumerate(x):
    xy.append((x, y[num]))
norm = 2
steps = ['ZI', 'IZ','ZI', 'IZ','ZI', 'IZ','ZI', 'IZ']
def get_ZI_line(xy_array=xy, Arai_steps=steps): # seems to work.  only possible problem is not excluding repeat ZI points
    ZI_points = []
    for num, point in enumerate(xy_array):
        if Arai_steps[num] == 'ZI':
            ZI_points.append(point)
    print "ZI points:", ZI_points
    ZI_line = 0.
    for num, point in enumerate(ZI_points[:-1]): # 
#        print "point", point
#        print "ZI_points[num + 1]", ZI_points[num+1]
        result = (numpy.sqrt( (point[0] - ZI_points[num + 1][0])**2 + (point[1] - ZI_points[num + 1][1])**2 ))
#        print "get sqrt of:", ( (point[0] - ZI_points[num + 1][0])**2 - (point[1] - ZI_points[num + 1][1])**2 )
#        print "(", point[0], "-", ZI_points[num + 1][0], ")**2  -  (", point[1], "-", ZI_points[num+1][1], ")**2"
#        print result
        ZI_line += result
    print ZI_line
    return { 'ZI_line': ZI_line, 'ZI_points': ZI_points }

# you need to work through this one by hand and figure it out, make sure it's right


xy_segment = ((1, 2),(2, 4),(5,6),(7,8),(8,10),(11,12),(13,14),(15,16),(17,18))
A_steps =     ['ZI',  'ZI', 'IZ',  'IZ', 'IZ',  'ZI',    'IZ',  'ZI',   'ZI']
x = [2, 8, 11, 13, 17]
y = [4, 10, 12, 14, 18]
real_ZI_line = get_ZI_line(xy_segment, A_steps)['ZI_line']
ref_triangles = [((2,4),(8,10),(11, 12)),((8,10),(11,12),(13,14)),((11,12),(13,14),(17,18))]
ref_midpoints = ['IZ', 'ZI', 'IZ']
# triangles and midpoints are both returned by get_triangles
def get_area_sum(triangles=ref_triangles, midpoints = ref_midpoints, ZI_line = real_ZI_line):  # this may work
    Area = 0
    for num, triangle in enumerate(triangles):
#        print triangle
        midpoint = midpoints[num]
        x_seg = [triangle[0][0], triangle[1][0], triangle[2][0]]
#        print "x_seg", x_seg  
        y_seg = [triangle[0][1], triangle[1][1], triangle[2][1]]
#        print "y_seg", y_seg
        lines = get_triangle_sides(x_seg, y_seg)  
#        print "lines", lines
        t = get_triangle(lines['L1'], lines['L2'], lines['L3'])
        t_area = t['triangle_A']
        sign = get_sign(triangle, midpoint)['sign']
#        print "sign", sign
        normed_triangle_area = (t_area * sign) / ZI_line
#        print "normed_triangle_area", normed_triangle_area
        Area += normed_triangle_area
#        print Area
    print "Area:", Area
    return Area


x_arai = [1, 2, 3, 4, 5, 6, 7]
y_arai = [1., .8, .6, .4, .2, .1, .05]
steps_arai = ['IZ', 'ZI', 'ZI', 'IZ', 'ZI', 'IZ', 'ZI']

def get_IZZI_MD(x_Arai=x_arai, y_Arai=y_arai, steps_Arai=steps_arai):
    norm = y_Arai[0] # initial NRM
   # def get_normed_points(point_array, norm): # good to go
    #return points
    
    x = get_normed_points(x_Arai, norm)
    y = get_normed_points(y_Arai, norm)

    #def get_xy_array(x_segment, y_segment): 
    #return xy_array

    xy_array = get_xy_array(x, y)

    #def get_triangles(xy_segment = xy_segment, Arai_steps = steps): # works!
    #return {'triangles': triangles, 'midpoints': midpoints}

    t = get_triangles(xy_array, steps_Arai)
    triangles = t['triangles']
    midpoints = t['midpoints']
    
    #def get_ZI_line(xy_array=xy, Arai_steps=steps): # seems to work.  only possible problem is not excluding repeat ZI points
    #return { 'ZI_line': ZI_line, 'ZI_points': ZI_points }

    ZI_line = get_ZI_line(xy_array, steps_Arai)['ZI_line']
    
   # def get_area_sum(triangles=ref_triangles, midpoints = ref_midpoints, ZI_line = real_ZI_line):  # this may work
    # calls get lines, get triangle, etc. for each triangle
    #return Area

    IZZI_MD = get_area_sum(triangles, midpoints, ZI_line)
    
    print "IZZI_MD:", IZZI_MD
    return IZZI_MD



