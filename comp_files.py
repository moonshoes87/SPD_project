#!/usr/bin/env python


# two files
# read each
# compare
# provide easy, readable diff

import sys

f1 = sys.argv.index('-f1')
f2 = sys.argv.index('-f2')
file1 = open(sys.argv[f1+1], 'rU')
file2 = open(sys.argv[f2+1], 'rU')

data1 = file1.readlines()[1]
data2 = file2.readlines()[1]

print data1[0], data1[1], data1[2]

l1 = data1.split('\t')
l2 = data2.split('\t')
print l1
print "-"
print l2


if True:
    for num, i in enumerate(l1[:-1]): # last item is ''.  this may change
        print "--"
        #print type(i), type(l2[num])
        #print i
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
        
        if n1 != n2:
            print i, "-----",  l2[num]
        else:
            print i.split()[0][:-1], " SAME!"
