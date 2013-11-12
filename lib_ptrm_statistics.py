#!/usr/bin/env python

import numpy



def get_n_ptrm(tmin, tmax, ptrm_temps, ptrm_starting_temps):
    """return number of ptrm_checks included in best fit segment.  excludes checks if temp exceeds tmax OR if starting temp exceeds tmax"""
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

def get_max_ptrm_check(ptrm_checks_segment, ptrm_checks, ptrm_x, t_Arai, x_Arai):
    diffs = []
    ptrm_checks = list(ptrm_checks)
    for check in ptrm_checks_segment:
        ptrm_ind = ptrm_checks.index(check)
        ptrm_check = ptrm_x[ptrm_ind]
        arai_ind = t_Arai.index(check)
        ptrm_orig = x_Arai[arai_ind]
        diff = ptrm_check / ptrm_orig
        diffs.append(diff)
    return max(diffs), max(diffs) * 100.

def get_delta_CK(max_ptrm_check, x_int):
    """Returns maximum difference produced by a ptrm check, normed by total TRM (x int of best fit line)"""
    return max_ptrm_check / x_int * 100.

def get_DRAT(delta_y_prime, delta_x_prime, max_ptrm_check):
    """Returns maximum difference produced by a ptrm check, normed by length of best_fit line"""
    L = numpy.sqrt(delta_x_prime**2 + delta_y_prime**2)
    DRAT = max_ptrm_check / L * 100
    return DRAT


    
