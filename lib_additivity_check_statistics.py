#!/usr/bin/env python

import numpy



def get_n_add(temps, starting_temps, tmax):
    incl_temps = []
    for num, temp in enumerate(temps):
        if temp <= tmax and starting_temps[num] <= tmax:
            incl_temps.append(temp)
    n_add = len(incl_temps)
    print incl_temps
    return incl_temps, n_add




def get_max_ptrm_check(ptrm_checks_included_temps, ptrm_checks_all_temps, ptrm_x, t_Arai, x_Arai):
    """sorts through included ptrm_checks and finds the largest ptrm check diff, the sum of the total diffs, and the percentage of the lar
gest check / original measurement at that temperature step"""
    diffs = []
    abs_diffs = []
    x_Arai_compare = []
    ptrm_compare = []
    check_percents = []
    ptrm_checks_all_temps = list(ptrm_checks_all_temps)
    for check in ptrm_checks_included_temps: # goes through each included temperature step
        ptrm_ind = ptrm_checks_all_temps.index(check) # indexes the number of the 
        ptrm_check = ptrm_x[ptrm_ind] # x value at that temperature step
        ptrm_compare.append(ptrm_check) #                                                                       
        arai_ind = t_Arai.index(check)
        ptrm_orig = x_Arai[arai_ind]
        x_Arai_compare.append(ptrm_orig)
        diff = ptrm_orig - ptrm_check
        diffs.append(diff)
        abs_diffs.append(abs(diff))
        check_percents.append((abs(diff) / ptrm_orig) * 100)
#    print "ptrm_checks_included_temps", ptrm_checks_included_
#    print "x_Arai_compare", x_Arai_compare
#    print "diffs", diffs                                                                                                
    max_diff = max(abs_diffs)
    check_percent = max(check_percents)
    sum_diffs = abs(sum(diffs))

