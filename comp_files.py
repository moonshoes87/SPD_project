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
#file1 = open(sys.argv[f1+1], 'rU')
#file2 = open(sys.argv[f2+1], 'rU')

#lines1 = file1.readlines()


def parse_file(file_path):
    """takes file and returns dictionary...."""
    file = open(file_path, 'rU')
    lines = file.readlines()
    data = []
    print len(lines)
    specimens = {}
    for line in lines[1:]: #
        d = line.split('\t')
        #print 'd        ', d
        data = []
        for i in d[:-1]: # empty space
            temp = i.split()
            data.append(temp[1])
        #print 'data', data
        specimens[data[0]] = data
    return specimens
    
specs1 = parse_file(path1)
specs2 = parse_file(path2)

print 'MINE:  ', specs1['ET1_279BS']
print 'GREIG\'S:   ', specs2['ET1_279BS']
print len(specs1['ET1_279BS'])
print len(specs2['ET1_279BS'])




if specs1.keys().sort() == specs2.keys().sort():
    print "yeya"
else:
    print "suck"

for specimen in specs1.keys():
    print "SPECIMEN: ", specimen
    data1 = specs1[specimen]
    print 'data1', data1
    data2 = specs2[specimen]
    for num, i in enumerate(data1):
        v1 = i
        v2 = data2[num]
        #print v1, v2
        str = False
        try:
            float(v1)
            r = round(float(v1), 1)
            n1 = r
            n2 = round(float(v2), 1)
            #print "ROUNDED", r
        except Exception as ex:
            #print ex
            str = True
        if not str:
            if isnan(float(v1)) and isnan(float(v2)):
             #   print " SAME!"
                break
        if str:
            if v1 == v2:
                pass
            else:
                print v1, "------", v2
              #  print "SAME"
        elif n1 != n2:
            print v1, "-----",  v2
        else:
            pass
            #print "same"
    print "--"
         
    

#print specs1#['ET1_283E']
#print '*********'
#print specs2#['ET1_283E']
