#!/usr/bin/env python


# two files
# read each
# compare
# provide easy, readable diff

import sys
from numpy import isnan

f1 = sys.argv.index('-f1')
f2 = sys.argv.index('-f2')
file1 = open(sys.argv[f1+1], 'rU')
file2 = open(sys.argv[f2+1], 'rU')

data1 = file1.readlines()[1]
data2 = file2.readlines()[1]


print file1
l1 = data1.split('\t')
print l1
print "-"
print file2
l2 = data2.split('\t')
print l2


if True:
    for num, i in enumerate(l1[:-1]): # last item is ''.  this may change
        print "--"
        #print type(i), type(l2[num])
        #print i
        str = False
        category = i.split()[0][:-1]
        n1 = i.split()[1]
        n2 = l2[num].split()[1]
        try:
            float(n1)
            r = round(float(n1), 1)
            n1 = r
            n2 = round(float(n2), 1)
            #print "ROUNDED", r
        except Exception as ex:
            print ex
            str = True
        if str == False:
            if isnan(n1) and isnan(n2):
                print category + " SAME!"
                break
        if n1 != n2:
            print i, "-----",  l2[num]
        else:
            print category +  " SAME!"
