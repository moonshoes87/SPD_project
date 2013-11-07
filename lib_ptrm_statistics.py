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
    return max(diffs) * 100.
#        return 0
    
