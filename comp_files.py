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
print data1[:300]
