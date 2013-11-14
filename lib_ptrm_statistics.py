#!/usr/bin/env python

import numpy



def get_n_ptrm(tmin, tmax, ptrm_temps, ptrm_starting_temps):
    """return number of ptrm_checks included in best fit segment.  excludes checks if temp exceeds tmax OR if starting temp exceeds tmax. also returns those ptrm_check temperatures"""
    # does not exclude ptrm checks that are less than tmin
    ptrm_checks_segment= []
    for num, check in enumerate(ptrm_temps):
        if check > tmax: #or check < tmin:
            pass
        elif ptrm_starting_temps[num] > tmax: # or ptrm_starting_temps[num] < tmin:
            pass
        else:
            ptrm_checks_segment.append(check)
    return len(ptrm_checks_segment), ptrm_checks_segment


def get_max_ptrm_check(ptrm_checks_included_temps, ptrm_checks_all_temps, ptrm_x, t_Arai, x_Arai): 
    """sorts through included ptrm_checks and finds the largest ptrm check diff, the sum of the total diffs, and the percentage of the largest check / original measurement at that temperature step"""
    diffs = []
    abs_diffs = []
    x_Arai_compare = []
    ptrm_compare = []
    check_percents = []
    ptrm_checks_all_temps = list(ptrm_checks_all_temps)
    for check in ptrm_checks_included_temps:
        ptrm_ind = ptrm_checks_all_temps.index(check)
        ptrm_check = ptrm_x[ptrm_ind]
        ptrm_compare.append(ptrm_check)
        arai_ind = t_Arai.index(check)
        ptrm_orig = x_Arai[arai_ind]
        x_Arai_compare.append(ptrm_orig)
        diff = ptrm_check - ptrm_orig
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


#DRATS = lib_ptrm.get_DRATS(self.ref_sum_ptrm_check, self.x_Arai, end)
    
