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

# fix it so that max_ptrm_check is right.  it should be subtraction, not division.  
def get_max_ptrm_check(ptrm_checks_included_temps, ptrm_checks_all_temps, ptrm_x, t_Arai, x_Arai): 
    diffs = []
    abs_diffs = []
    x_Arai_compare = []
    ptrm_compare = []
#    print "ptrm_checks_included_temps", ptrm_checks_included_temps
    ptrm_checks_all_temps = list(ptrm_checks_all_temps)
#    print "ptrm_x", ptrm_x
#    print "ptrm_checks_all_temps", ptrm_checks_all_temps
    for check in ptrm_checks_included_temps:
#        print "check", check
        ptrm_ind = ptrm_checks_all_temps.index(check)
#        print "ptrm_ind", ptrm_ind
        ptrm_check = ptrm_x[ptrm_ind]
        ptrm_compare.append(ptrm_check)
        arai_ind = t_Arai.index(check)
        ptrm_orig = x_Arai[arai_ind]
        x_Arai_compare.append(ptrm_orig)
#        print "ptrm_check", ptrm_check, "ptrm_orig", ptrm_orig
        diff = ptrm_check - ptrm_orig
        diffs.append(diff)
        abs_diffs.append(abs(diff))
#    print "diffs", diffs
#    print "ptrm_checks", ptrm_compare
#    print "x_Arai_compare", x_Arai_compare
    print "abs_diffs", abs_diffs
    print "diffs", diffs
    max_diff = max(abs_diffs)
    diff_ind = abs_diffs.index(max_diff)
    print diff_ind
    norm_x = x_Arai_compare[diff_ind]
    print x_Arai_compare, "x_Arai_compare"
    print norm_x
    check_percent = max_diff / norm_x * 100
    return max_diff, sum(diffs), check_percent

def get_check_percent(max_ptrm_diff, x_at_that_step):
    pass


def get_delta_CK(max_ptrm_check, x_int):
    """Returns maximum difference produced by a ptrm check, normed by total TRM (x int of best fit line)"""
    return max_ptrm_check / x_int * 100.

def get_DRAT(delta_y_prime, delta_x_prime, max_ptrm_check):
    """Returns maximum difference produced by a ptrm check, normed by length of best_fit line"""
    L = numpy.sqrt(delta_x_prime**2 + delta_y_prime**2)
    DRAT = max_ptrm_check / L * 100
    return DRAT


    
