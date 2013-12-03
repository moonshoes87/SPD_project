#!/usr/bin/env python

import numpy
import lib_directional_statistics as lib_direct

numpy.set_printoptions(precision=15)


def get_n_ptrm(tmin, tmax, ptrm_temps, ptrm_starting_temps):
    """return number of ptrm_checks included in best fit segment.  excludes checks if temp exceeds tmax OR if starting temp exceeds tmax. also returns those ptrm_check temperatures"""
    # does not exclude ptrm checks that are less than tmin
    ptrm_checks_included_temps= []
    for num, check in enumerate(ptrm_temps):
        if check > tmax: 
            pass
        elif ptrm_starting_temps[num] > tmax: # or ptrm_starting_temps[num] < tmin:
            pass
        else:
            ptrm_checks_included_temps.append(check)
    return len(ptrm_checks_included_temps), ptrm_checks_included_temps


def get_max_ptrm_check(ptrm_checks_included_temps, ptrm_checks_all_temps, ptrm_x, t_Arai, x_Arai): 
    """sorts through included ptrm_checks and finds the largest ptrm check diff, the sum of the total diffs, and the percentage of the largest check / original measurement at that temperature step"""
    diffs = []
    abs_diffs = []
    x_Arai_compare = []
    ptrm_compare = []
    check_percents = []
    ptrm_checks_all_temps = list(ptrm_checks_all_temps)
    for check in ptrm_checks_included_temps: # goes through each included temperature step
        ptrm_ind = ptrm_checks_all_temps.index(check) # indexes the number of the check
        ptrm_check = ptrm_x[ptrm_ind] # x value at that temperature step
        ptrm_compare.append(ptrm_check) # 
        arai_ind = t_Arai.index(check) 
        ptrm_orig = x_Arai[arai_ind]
        x_Arai_compare.append(ptrm_orig)
        diff = ptrm_orig - ptrm_check
        diffs.append(diff)
        abs_diffs.append(abs(diff))
        check_percents.append((abs(diff) / ptrm_orig) * 100)
#    print "ptrm_checks_included_temps", ptrm_checks_included_temps
#    print "ptrm_compare", ptrm_compare
#    print "x_Arai_compare", x_Arai_compare
#    print "diffs", diffs
    max_diff = max(abs_diffs)
    check_percent = max(check_percents)
    sum_diffs = abs(sum(diffs))
    sum_abs_diffs = sum(abs_diffs)
    return max_diff, sum_diffs, check_percent, sum_abs_diffs # gives the largest difference between an original trm measurement and a ptrm check, the sum of those differences, the percentage diff of the largest difference, the absolute value sum of the diffs

def get_delta_CK(max_ptrm_check, x_int):
    """Returns maximum difference produced by a ptrm check, normed by total TRM (x int of best fit line)"""
    return (max_ptrm_check / x_int) * 100.

def get_DRAT(delta_x_prime, delta_y_prime, max_ptrm_check):
    """Returns maximum difference produced by a ptrm check, normed by length of best_fit line"""
    L = numpy.sqrt(delta_x_prime**2 + delta_y_prime**2)
    DRAT = (max_ptrm_check / L) * 100
    return DRAT, L

def get_max_DEV(delta_x_prime, max_ptrm_check):
    return (max_ptrm_check / delta_x_prime) * 100.

def get_CDRAT(L, sum_ptrm_checks, sum_abs_ptrm_checks):
    CDRAT = (sum_ptrm_checks / L) * 100.
    CDRAT_prime = (sum_abs_ptrm_checks / L) * 100.
    return CDRAT, CDRAT_prime

def get_DRATS(sum_ptrm_checks, sum_abs_ptrm_checks, x_Arai, end):
    DRATS = (sum_ptrm_checks / x_Arai[end]) * 100.
    DRATS_prime = (sum_abs_ptrm_checks / x_Arai[end]) * 100.
    print "x_Arai[end]", x_Arai[end]
    return DRATS, DRATS_prime

def get_mean_DRAT(sum_ptrm_checks, sum_abs_ptrm_checks, n_pTRM, L):
    mean_DRAT = ((1. / n_pTRM) * (sum_ptrm_checks / L)) * 100
    mean_DRAT_prime = ((1./ n_pTRM) * (sum_abs_ptrm_checks / L)) * 100
    return mean_DRAT, mean_DRAT_prime

def get_mean_DEV(sum_ptrm_checks, sum_abs_ptrm_checks, n_pTRM, delta_x_prime):
    mean_DEV = ((1. / n_pTRM) * (sum_ptrm_checks / delta_x_prime)) * 100
    mean_DEV_prime= ((1. / n_pTRM) * (sum_abs_ptrm_checks / delta_x_prime)) * 100
    return mean_DEV, mean_DEV_prime
#    mean_DEV, mean_DEV_prime = lib_ptrm.get_mean_DEV(self.ref_sum_ptrm_check, self.ref_sum_abs_ptrm_check, self.ref_n, self.delta_x_prime)

#DRATS = lib_ptrm.get_DRATS(self.ref_sum_ptrm_check, self.x_Arai, end)
    
def get_delta_pal_vectors(PTRMS, PTRM_Checks):
    """ takes in PTRM data in this format: [temp, dec, inc, moment, ZI or IZ] -- and PTRM_check data in this format: [temp, dec, inc, moment].  Returns them in vector form (cartesian). """
    if type(PTRMS) != numpy.ndarray:
        PTRMS = numpy.array(PTRMS)
    if type(PTRM_Checks != numpy.ndarray):
        PTRM_Checks = numpy.array(PTRM_Checks)
    TRM_1 = lib_direct.dir2cart(PTRMS[0,1:3])
    PTRMS_cart = lib_direct.dir2cart(PTRMS[1:,1:3])
    checks = lib_direct.dir2cart(PTRM_Checks[:,1:3])
    return PTRMS_cart, checks, TRM_1

def get_diffs(ptrms_vectors, ptrm_checks_vectors, ptrms_orig, checks_orig):  
    """presumes ptrms_orig & checks orig have format [[temp, dec, inc, moment, zi/iz], ...].  ptrms_vectors and ptrm_checks_vectors are cartesian [[x,y,z],...].  requires both vector and original format of ptrms and checks for correct temperature indexing.  takes these in and returns diffs and C"""
#    print "ptrms_vectors", ptrms_vectors
#    print "ptrm_checks_vectors", ptrm_checks_vectors
#    print "ptrms original", ptrms_orig
#    print "checks original", checks_orig
    if len(ptrms_vectors) == len(ptrm_checks_vectors):
        diffs = ptrms_vectors - ptrm_checks_vectors
    else:
        diffs = numpy.zeros((len(ptrms_vectors), 3))
        check_num = 0
        for num, ptrm in enumerate(ptrms_orig):
            if num == 0:
#                print "first is zero"
                diff = [0,0,0]
            elif len(checks_orig) <= check_num:  # if there are ptrms at higher temps than checks, this breaks the cycle
                break
            elif ptrms_orig[num][0] == checks_orig[check_num][0]:
#                print "ptrm", ptrms_orig[num], "check", checks_orig[check_num]
                diff = ptrms_vectors[num-1] - ptrm_checks_vectors[check_num]
                check_num += 1
            else:
#                print "ptrm", ptrms_orig[num], "check", checks_orig[check_num]
                diff = [0,0,0]
            diffs[num-1] = diff
#            print "diff:", diff
    C = numpy.cumsum(diffs, axis=0)
    return diffs, C


def get_TRM_star(PTRMS_cart, C, TRM_1):
    TRM1 = TRM_1.reshape((1,3)) # ensures that TRM_1 is compatible to be concatenated
    TRMS_adj = PTRMS_cart + C
    TRM_star = numpy.concatenate((TRM1, TRMS_adj))
    x_star = numpy.zeros(len(TRM_star))
    for num, trm in enumerate(TRM_star):
        x_star[num] = numpy.linalg.norm(trm)
    return TRM_star, x_star

def get_b_star(x_star, y_err, y_mean):
    """get corrected x segment and x_mean"""
    x_star_mean = numpy.mean(x_star)
    x_err = x_star - x_star_mean
    b_star = -1* numpy.sqrt( sum(y_err**2) / sum(x_err**2) )  # averaged slope 
    return b_star

def get_delta_pal(b, b_star):
    delta_pal = numpy.abs((b - b_star) / b) * 100
    return delta_pal

def get_full_delta_pal(PTRMS, PTRM_Checks, y_err, y_mean, b):
#    return 0
    PTRMS_cart, checks, TRM_1 = get_delta_pal_vectors(PTRMS, PTRM_Checks)
#    print "PTRMS_Cart", PTRMS_cart
    diffs, C = get_diffs(PTRMS_cart, checks, PTRMS, PTRM_Checks)
#    print "C", C
    TRM_star, x_star = get_TRM_star(PTRMS_cart, C, TRM_1)
#    print "x_star", x_star
#    print type(x_star)
    b_star = get_b_star(x_star, y_err, y_mean)
    delta_pal = get_delta_pal(b, b_star)
    return delta_pal


def get_segments(ptrms, ptrm_checks, tmax):
    ptrms_included = []
    checks_included = []
    ptrms = numpy.array(ptrms)
    for ptrm in ptrms:
        if ptrm[0] <= tmax:
            ptrms_included.append(ptrm)
    for check in ptrm_checks:
        if check[0] <= tmax:
            checks_included.append(check)
    print "checks", ptrm_checks
    print "checks_included", checks_included
    return ptrms_included, checks_included

    
#    def get_delta_pal_vectors(PTRMS, PTRM_Checks):
    #return PTRMS_cart, checks, TRM_1

    #def get_diffs(ptrms_vectors, ptrm_checks_vectors, ptrms_orig, checks_orig):  
#    return diffs, C

#   def get_TRM_star(PTRMS_cart, C, TRM_1):
#    return TRM_star, x_star

#def get_b_star(x_star, y_err, y_mean):
    #return b_star

#def get_delta_pal(b, b_star):
    # return delta_pal

            

# york b code
#    x_err = x_segment - x_mean
#    y_err = y_segment - y_mean
#    york_b = -1* sqrt( sum(y_err**2) / sum(x_err**2) )  # averaged slope 



