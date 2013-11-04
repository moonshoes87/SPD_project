#!/usr/bin/env python 

import numpy

    
"""
function Par = TaubinSVD(XY)

%--------------------------------------------------------------------------
%  
%     Algebraic circle fit by Taubin
%      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
%                  Space Curves Defined By Implicit Equations, With 
%                  Applications To Edge And Range Image Segmentation",
%      IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
    XY = [[1, 2], [2,3], [4, 5]]
    x_i = XY[i][0]
    y_i = XY[i][1]
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this is a version optimized for stability, not for speed
%
%--------------------------------------------------------------------------

"""

xy = [[0.,-2.], [2.5,3.], [4., 5.]]

def TaubinSVD(XY = numpy.array(xy)):
    print "XY", XY
    X = XY[:,0] - numpy.mean(XY[:,0]) # norming points by x avg
    Y = XY[:,1] - numpy.mean(XY[:,1]) # norming points by y avg
    centroid = [numpy.mean(XY[:,0]), numpy.mean(XY[:,1])]
    print "centroid", centroid
#    XY[:,0] = XY[:, 0] - centroid[0] # from each x, subtract x_avg
#    XY[:,1] = XY[:, 1] - centroid[1] # from each y, subtract y_avg
    print "XY", XY
    print "X", X, "Y", Y
    # Z is correct
    Z = X * X + Y * Y  # Z = X.*X + Y.*Y; .*   # in matlab, .* is equivalent to *, and * is equivalent to numpy.dot
    print "Z", Z
    Zmean = numpy.mean(Z)
    print "Zmean", Zmean
    Z0 = (Z - Zmean) / (2. * numpy.sqrt(Zmean))
    ZXY = numpy.array([Z0, X, Y])
    print "ZXY", ZXY
    U, S, V = numpy.linalg.svd(ZXY, full_matrices=False) # svd(X, 0) in original documentation.  however, belive that full_matrices=False accomplishes same
    print "U", U
    print "S", S
    print "V", V
    A = V[:,2]
    print "A", A
    A[0] = A[0] / (2. * numpy.sqrt(Zmean))
    print "adjusted A", A
    print "A", A, "other part", (-1. * Zmean * A[0])
    A = numpy.concatenate([A, [(-1. * Zmean * A[0])]], axis=0)
    print "final A", A
    #          -(A(2:3))'/A(1)/2+centroid
    a, b = (-1 * A[1:2]) / A[0] / 2 + centroid # but should be transposed, somewhere???.  syntax is problematic.  
    print "a,b", a,b
    #         sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2];
    r = numpy.sqrt(A[1]*A[1]+A[2]*A[2]-4*A[0]*A[3])/abs(A[0])/2;
#    return A, V
    return a,b,r, XY[:,0], XY[:,1]
    
#    return (-1 * numpy.transpose((A[1:2]))) / A[0] / 2 + centroid
    
ignore = """
#Transpose it to a column vector -- use ' 
centroid = mean(XY);   % the centroid of the data set
X = XY(:,1) - centroid(1);  %  centering data
Y = XY(:,2) - centroid(2);  %  centering data
Z = X.*X + Y.*Y;
Zmean = mean(Z);
Z0 = (Z-Zmean)/(2*sqrt(Zmean));
ZXY = [Z0 X Y];
[U,S,V]=svd(ZXY,0);
A = V(:,3);
A(1) = A(1)/(2*sqrt(Zmean));
A = [A ; -Zmean*A(1)];  # the semi-colon means concatenate((a,b), axis=0)
Par = [-(A(2:3))'/A(1)/2+centroid , sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2];

end  
"""
