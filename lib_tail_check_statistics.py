#!/usr/bin/env python                                                                                            
import numpy


def get_n_tail(tmax, tail_temps):
    """determines number of included tail checks in best fit segment"""
    t_index = 0
    adj_tmax = 0
    try:
        t_index = list(tail_temps).index(tmax)
    except: # finds correct tmax if there was no tail check performed at tmax
        for temp in tail_temps:
            if temp <= tmax:
                adj_tmax = temp
        t_index = list(tail_temps).index(adj_tmax)
    incl_temps = tail_temps[0:t_index+1] # b/c not inclusive
    return len(incl_temps) #, incl_temps


def get_max_tail_check(y_Arai, y_tail, t_Arai, tail_temps, n_tail):
#tail_check_max, tail_check_diffs = lib_tail.get_tail_check_max(self.y_Arai, self.y_tail, self.t_Arai, self.tail_temps, self.ref_n_tail)
    tail_compare = []
    y_Arai_compare = []
    print "n_tail", n_tail
    for temp in tail_temps[:n_tail]:
        print "temp", temp
        tail_index = tail_temps.index(temp)
        tail_check = y_tail[tail_index]
        tail_compare.append(tail_check)
        arai_index = t_Arai.index(temp)
        print "arai_index", arai_index
        nrm_orig = y_Arai[arai_index]
        y_Arai_compare.append(nrm_orig)
    print "tail compare", tail_compare
    print "y_Arai_compare", y_Arai_compare
    diffs = numpy.array(y_Arai_compare) - numpy.array(tail_compare)
    abs_diffs = abs(diffs)
    max_check = max(abs_diffs)
    print "diffs", diffs
    return max_check, diffs


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
    max_diff = max(abs_diffs)
    check_percent = max(check_percents)
    sum_diffs = abs(sum(diffs))
    sum_abs_diffs = sum(abs_diffs)
