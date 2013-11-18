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
        tail_index = list(tail_temps).index(temp)
        tail_check = y_tail[tail_index]
        tail_compare.append(tail_check)
        arai_index = list(t_Arai).index(temp)
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

