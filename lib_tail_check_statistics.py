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

