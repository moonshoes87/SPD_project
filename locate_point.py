#! /usr/bin/env/python

#line 1:
slope1 = -1
y_int1 = 5
x_int1 = 5

def return_y(slope, y_int, x):
    """
    slope, y_int, x
    """
    y = (slope * x) + y_int
    return y


slope2 = -1
y_int2 = 1.5
x_int2 = 1.5

x_max = x_int1
y_max = y_int1


def locate_point(x, y):
    fail = False
    if x > x_max:
        print "x > x_max"
    if y > y_max:
        print "y > y_max"
    if x < 0:
        print "x < 0"
    if y < 0:
        print "y < 0 "
    biggest_y = return_y(slope1, y_int1, x)
    print "biggest_y", biggest_y
    smallest_y = return_y(slope2, y_int2, x)
    print "smallest_y", smallest_y
    if y > biggest_y or y < smallest_y:
        print "y > biggest_y or y < smallest_y:"



