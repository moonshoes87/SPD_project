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
    for num, i in enumerate(l1):
        print "--"
        if i != l2[num]:
            print i, "-----",  l2[num]
