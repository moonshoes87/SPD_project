#!/usr/bin/env python

import numpy



def get_n_ptrm(tmin, tmax, ptrm_temps, ptrm_starting_temps):
    """return number of ptrm_checks included in best fit segment.  excludes checks if temp exceeds tmax OR if starting temp exceeds tmax"""
    ptrm_steps = []
    for num, check in enumerate(ptrm_temps):
        if check > tmax or check < tmin:
            pass
        elif ptrm_starting_temps[num] > tmax or ptrm_starting_temps[num] < tmin:
            pass
        else:
            ptrm_steps.append(check)
    return len(ptrm_steps)

#    for num, temp in enumerate(ptrm_checks_temperatures): # seems 
#        if temp >= tmin and temp <= tmax: # if temp is within selected range  
#            if ptrm_checks_starting_temperatures[num] >= tmin and ptrm_checks_starting_temperatures[num] <= tmax: # and also if it was not done after an out-of-range temp
#                x = x_ptrm_check[num]
#                y = y_ptrm_check[num]
#                points.append((x, y))
