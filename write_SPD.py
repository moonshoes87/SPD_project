#!/usr/bin/env python


import numpy
import spd



# maybe automate a way to go through multiple magic_measurements.txt files in a directory
# or just have as input all the files and make a "gui" object with each
#  get list of specimens.  go through them with tmin = min possible and tmax = max possible
# use list(s) to get order of tabbed file nice
# make sure you use the correct files (the ones that you have comparisons for)


in_files = ['magic_measurements_ET1_Biggin_2007.txt']
outfile = 'ET1_Biggin_2007_output.txt'

        
basic_stats = ['s', 'specimen_n', 'tmin', 'tmax'] # start, end  # removed these two because no equivalent in Greig's code

arai_plot_stats = ['specimen_b', 'specimen_b_sigma', 'B_lab', 'B_anc', 'B_anc_sigma', 'specimen_YT', 'specimen_XT', 'specimen_vds',  'delta_x_prime', 'delta_y_prime', 'specimen_f', 'specimen_fvds', 'FRAC', 'specimen_b_beta', 'specimen_g', 'GAP-MAX', 'specimen_q', 'specimen_w', 'specimen_k', 'SSE', 'SCAT', 'R_corr2', 'R_det2', 'Z', 'Zstar', 'IZZI_MD']  #  'x_prime', 'y_prime',  # removed these two because they are listy and it sucks

directional_stats = ['Dec_Free', 'Dec_Anc', 'Inc_Free', 'Inc_Anc', 'MAD_Free', 'MAD_Anc', 'alpha', 'theta', 'DANG', 'NRM_dev', 'gamma']

ptrm_stats = ['n_ptrm', 'max_ptrm_check_percent', 'delta_CK', 'DRAT', 'max_DEV', 'CDRAT', 'CDRAT_prime', 'DRATS', 'DRATS_prime', 'mean_DRAT', 'mean_DRAT_prime', 'mean_DEV', 'mean_DEV_prime', 'delta_pal']

tail_stats = ['n_tail', 'DRAT_tail', 'delta_TR', 'MD_VDS']

additivity_stats = ['n_add', 'delta_AC']

long_list = arai_plot_stats + directional_stats + ptrm_stats + tail_stats + additivity_stats




import new_lj_thellier_gui_spd as tgs

print "starting thingee"


# create tab file
out = open(outfile, 'w')
for stat in basic_stats:
    out.write(stat + '\t')
for stat in long_list:
    out.write(stat + '\t')
out.write('\n')


# iterate through list of magic_measurements files and fill outfile

for f in in_files:
    gui = tgs.Arai_GUI(f)
    data = gui.Data
    specimen_names = gui.Data.keys()
    for s in specimen_names:
#        spec = spd.PintPars(gui.Data, s, gui.Data[s]['t_Arai'][0], gui.Data[s]['t_Arai'][6])
        spec = spd.PintPars(gui.Data, s, gui.Data[s]['t_Arai'][1], gui.Data[s]['t_Arai'][3])
        print spec.s, spec.tmin, spec.tmax
        spec.calculate_all_statistics()
        out.write("s: {} \t n: {} \t Tmin: {} \t Tmax: {} \t".format(str(spec.s), str(spec.pars['specimen_n']), str(spec.tmin_K), str(spec.tmax_K)))
        for stat in long_list:
            if type(spec.pars[stat]) == numpy.ndarray:  # catches arrays to prevent extra newlines being auto-inserted
                out.write(str(stat) + ": " + numpy.array_str(spec.pars[stat], max_line_width=10000000) + '\t')
            else:
                out.write(str(stat) + ": " + str(spec.pars[stat]) + '\t')
        out.write('\n')
out.close()

