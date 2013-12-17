#!/usr/bin/env python



import spd



# maybe automate a way to go through multiple magic_measurements.txt files in a directory
# or just have as input all the files and make a "gui" object with each
#  get list of specimens.  go through them with tmin = min possible and tmax = max possible
# use list(s) to get order of tabbed file nice
# make sure you use the correct files (the ones that you have comparisons for)


in_files = ['magic_measurements.txt']
outfile = 'zee_output.txt'

        
basic_stats = ['s', 'specimen_n', 'start', 'end', 'tmin', 'tmax']

arai_plot_stats = ['specimen_n', 'start', 'end', 'tmin', 'tmax', 'specimen_b', 'specimen_b_sigma', 'B_anc', 'B_anc_sigma', 'specimen_YT', 'specimen_XT', 'specimen_vds', 'x_prime', 'y_prime', 'delta_x_prime', 'delta_y_prime', 'specimen_f', 'specimen_fvds', 'FRAC', 'specimen_b_beta', 'specimen_g', 'GAP-MAX', 'specimen_q', 'specimen_w', 'specimen_k', 'SSE', 'SCAT', 'R_corr2', 'R_det2', 'Z', 'Zstar', 'IZZI_MD']

directional_stats = ['Dec_free', 'Dec_anc', 'Inc_Free', 'Inc_Anc', 'MAD_Free', 'MAD_Anc', 'alpha', 'theta', 'DANG', 'NRM_dev', 'gamma']

ptrm_stats = ['n_ptrm', 'max_ptrm_check_percent', 'delta_CK', 'DRAT', 'max_DEV', 'CDRAT', 'DRATS', 'DRATS_prime', 'mean_DRAT', 'mean_DRAT_prime', 'mean_DEV', 'mean_DEV_prime', 'delta_pal']

tail_stats = ['n_tail', 'DRAT_tail', 'delta_TR', 'MD_VDS']

additivity_stats = ['n_add', 'Delta_AC']

long_list = basic_stats + arai_plot_stats + directional_stats + ptrm_stats + tail_stats + additivity_stats




import new_lj_thellier_gui_spd as tgs

print "starting thingee"

for f in in_files:
    gui = tgs.Arai_GUI(f)
    out = open(outfile, 'a')
    specimens = gui.Data.keys()
    for s in specimens:
        spec = gui.Data[s]
        for stat in arai_plot_stats:
            out.write(str(spec.pars[stat]) + '\t')
        out.write('\n')
