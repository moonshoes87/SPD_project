#!/usr/bin/env python


# two files
# read each
# compare
# provide easy, readable diff

import sys
from numpy import isnan



f1 = sys.argv.index('-f1')
f2 = sys.argv.index('-f2')
path1 = sys.argv[f1+1]
path2 = sys.argv[f2+1]


def parse_file(file_path):
    """takes file and returns dictionary...."""
    file = open(file_path, 'rU')
    lines = file.readlines()
    categories = []
    data = []
    specimens = {}
    for line in lines[1:]: #
        d = line.split('\t')
        #print 'd        ', d
        data = []
        for i in d[:-1]: # empty space
            temp = i.split()
            categories.append(temp[0])
            data.append(temp[1])
        specimens[data[0]] = data
    return specimens, categories
    

#print 'MINE:  ', specs1['ET1_279BS']
#print 'GREIG\'S:   ', specs2['ET1_279BS']
#print len(specs1['ET1_279BS'])
#print len(specs2['ET1_279BS'])



def compare_all(categories, specs1, specs2):

    if specs1.keys().sort() == specs2.keys().sort():
        pass
    else:
        raise NameError('different specimens detected')

    for specimen in specs1.keys():
        print "SPECIMEN: ", specimen
        data1 = specs1[specimen]
        data2 = specs2[specimen]
        for num, i in enumerate(data1):
            v1 = i
            v2 = data2[num]
            #print categories1[num]
            #print v1, v2
            str = False
            try:
                float(v1)
                n1 = round(float(v1), 1)
                n2 = round(float(v2), 1)
                #print "ROUNDED", r
            except ValueError as ex: # could not convert string to float
                #print ex
                str = True
            if not str:
                if isnan(float(v1)) and isnan(float(v2)):
                    #   print " SAME!"
                    break
            if str:
                if categories1[num] == 'SCAT:':
                    if bool(v1) == bool(v2):
                        pass
                    else:
                        print categories1[num]
                        print v1, "------", v2
                elif v1 == v2: # turn this back to pass!
                    print categories1[num]
                    print v1, "-----",  v2
                    #pass
                else:
                    print categories1[num]
                    print v1, "------", v2
                    #  print "SAME"
            elif n1 != n2:
                print categories1[num]
                print v1, "-----",  v2
            else:
                pass
                #print "same"
        print "--"



def compare_one(stat, categories, specs1, specs2):
    print "comparing one"
    print stat
    i = categories.index(stat + ':')
    key = specs1.keys()[-1]
    print specs1[key][i]
    print specs2[key][i]
    return True


specs1, categories1 = parse_file(path1)
specs2, categories2 = parse_file(path2)
try:
    single = sys.argv.index('-s')
    stat = sys.argv[single+1]
    compare_one(stat, categories1, specs1, specs2)
except ValueError as ex:
    print ex
    compare_all(categories1, specs1, specs2)
         
    
lori_only = ['y_prime', 'x_prime', 'delta_y_prime', 'delta_x_prime']
greig_only = ['Line_Len']
#print lori_only
#print greig_only
