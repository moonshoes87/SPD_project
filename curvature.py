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
def TaubinSVD(XY):
    centroid = numpy.mean(XY)
    XY[:, 0] -= centroid[0] # from each x, subtract x_avg
    XY[:, 1] -= centroid[1] # from each y, subtract y_avg
    X = XY[:,0]
    Y = XY[:,1]
    Z = X * X + Y * Y  # Z = X.*X + Y.*Y; .*   # in matlab, .* is equivalent to *, and * is equivalent to numpy.dot
    Zmean = numpy.mean(Z)
    Z0 = (Z - Zmean) / (2. * numpy.sqrt(Zmean))
    ZXY = numpy.array([Z0, X, Y])
    U, S, V = numpy.linalg.svd(ZXY, full_matrices=False) # svd(X, 0) in original documentation.  however, belive that full_matrices=False accomplishes same
    A = V[:,2]
    A[0] = A[0] / (2. * numpy.sqrt(Zmean))
    A = [[A], [-1. * Zmean * A[0]]]
    return 
    


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
A = [A ; -Zmean*A(1)];
Par = [-(A(2:3))'/A(1)/2+centroid , sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2];

end  
