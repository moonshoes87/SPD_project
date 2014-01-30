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
    """
    input: y_Arai, y_tail, t_Arai, tail_temps, n_tail
    output: max_check, diffs
    """
    tail_compare = []
    y_Arai_compare = []
    for temp in tail_temps[:n_tail]:
        tail_index = list(tail_temps).index(temp)
        tail_check = y_tail[tail_index]
        tail_compare.append(tail_check)
        arai_index = list(t_Arai).index(temp)
        nrm_orig = y_Arai[arai_index]
        y_Arai_compare.append(nrm_orig)
    diffs = numpy.array(y_Arai_compare) - numpy.array(tail_compare)
    abs_diffs = abs(diffs)
    max_check = max(abs_diffs)
    return max_check, diffs

def get_DRAT_tail(max_check, L):
    """
    input: tail_check_max, best fit line length
    output: DRAT_tail
    """
    DRAT_tail = (max_check / L) * 100.
    return DRAT_tail

def get_delta_TR(tail_check_max, y_int):
    """
    input: tail_check_max, y_intercept
    output: delta_TR
    """
    delta_TR = (tail_check_max / y_int) * 100.
    return delta_TR

def get_MD_VDS(tail_check_max, vds):
    """
    input: tail_check_max, vector difference sum
    output: MD_VDS
    """    
    MD_VDS = (tail_check_max / vds) * 100
    return MD_VDS
