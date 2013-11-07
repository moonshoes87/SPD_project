#!/usr/bin/env python 

import numpy
#import spd
#import lib_arai_plot_statistics as lib_arai

#thing = spd.thing
#  def get_xy_array(x_segment, y_segment):
#xy_real = lib_arai.get_xy_array(thing.x_Arai_segment, thing.y_Arai_segment)
#xy_real = [[0.12182379795523131, 0.98795180722891562], [0.1876743181924827, 0.95783132530120485], [0.21711848591334654, 0.96987951807228912], [0.32091411888460553, 0.98192771084337338], [0.48291503240724426, 0.9337349397590361], [0.72423703304017273, 0.90361445783132521], [1.0313987625711432, 0.81325301204819278]]    
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

xy_real = [[0.12182379795523131, 0.98795180722891562], [0.1876743181924827, 0.95783132530120485], [0.21711848591334654, 0.96987951807228912], [0.32091411888460553, 0.98192771084337338], [0.48291503240724426, 0.9337349397590361], [0.72423703304017273, 0.90361445783132521], [1.0313987625711432, 0.81325301204819278]]    

# coming out fairly different from get_curvature in spd.  need matlab

xy = [[0.,-2.], [2.5,3.], [4., 5.]]

def TaubinSVD(XY = xy):
    XY = numpy.array(XY)
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
    print "Z0", Z0
    ZXY = numpy.array([Z0, X, Y]).T
    print "ZXY", ZXY
    U, S, V = numpy.linalg.svd(ZXY, full_matrices=False) # svd(X, 0) in original documentation.  however, belive that full_matrices=False accomplishes same
    print "U", U
    print "S", S
    print "V", V
    V = V.transpose()
    A = V[:,2]
    print "A", A
    A[0] = A[0] / (2. * numpy.sqrt(Zmean))
    print "adjusted A", A
    print "A", A, "other part", (-1. * Zmean * A[0])
    A = numpy.concatenate([A, [(-1. * Zmean * A[0])]], axis=0)
    print "final A", A
    #          -(A(2:3))'/A(1)/2+centroid
    a, b = (-1 * A[1:3]) / A[0] / 2 + centroid # but should be transposed, somewhere???.  syntax is problematic.  
    # A[1:3].conj().transpose()
    print "a,b", a,b
    #         sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2];
    r = numpy.sqrt(A[1]*A[1]+A[2]*A[2]-4*A[0]*A[3])/abs(A[0])/2;
#    return A, V
    return { 'a':a,'b': b, 'r': r } #, XY[:,0], XY[:,1]
    
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



xy = numpy.array([[1,3],[0,2], [3,7], [8,10], [9,12]])
par = [1,2,7]

def VarCircle(XY = xy, Par = par):  # must have at least 4 sets of xy points or else division by zero occurs
    """
    computing the sample variance of distances from data points (XY) to the circle Par = [a b R]
    """
    print "varcircle"
    print "XY", XY
    n = len(XY)
    print "n", n
    Dx = XY[:,0] - Par[0]
    Dy = XY[:,1] - Par[1]
    D = numpy.sqrt(Dx * Dx + Dy * Dy) - Par[2]
    result = numpy.dot(D, D)/(n-3)
    print n, Dx, Dy
    print D
    print result
    print "done varcircle"
    return result


xy = numpy.array([[1,1],[0.,-2.], [.5,3.], [4., 5.]])
par_ini = [21.2500,   -9.5000,   22.5347]

def LMA(XY=xy,ParIni=par_ini):
    """
    %     Geometric circle fit (minimizing orthogonal distances)  
    %     based on the Levenberg-Marquardt scheme in the
    %     "algebraic parameters" A,B,C,D  with constraint B*B+C*C-4*A*D=1
    %        N. Chernov and C. Lesort, "Least squares fitting of circles",
    %        J. Math. Imag. Vision, Vol. 23, 239-251 (2005)
    %     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)             
    %             ParIni = [a b R] is the initial guess (supplied by user)
    %     Output: Par = [a b R] is the fitting circle:                                               
    %                           center (a,b) and radius R                                                                  
    %                                                                                                     
    """
    factorUp=10
    factorDown=0.04
    lambda0=0.01
    epsilon=0.000001
    IterMAX = 50
    AdjustMax = 20
    Xshift=0  
    Yshift=0  
    dX=1  
    dY=0;                                                                                    
    n = len(XY);      # number of data points  # checked
#    print n

    anew = ParIni[0] + Xshift
    bnew = ParIni[1] + Yshift
    Anew = 1./(2.*ParIni[2])                                                                              
    aabb = anew*anew + bnew*bnew    
    Fnew = (aabb - ParIni[2]*ParIni[2])*Anew # checked
#    print "Fnew", Fnew
    Tnew = numpy.arccos(-anew/numpy.sqrt(aabb)) # checked
#    print "Tnew", Tnew

    if bnew > 0:
        Tnew = 2*numpy.pi - Tnew
#    print "XY & ParIni", XY, ParIni
    VarNew = VarCircle(XY,ParIni) # checked
#    print VarNew 




ignore_me_too = """


function Par = LMA(XY,ParIni)
%--------------------------------------------------------------------------                                                %     Geometric circle fit (minimizing orthogonal distances)  
%     based on the Levenberg-Marquardt scheme in the
%     "algebraic parameters" A,B,C,D  with constraint B*B+C*C-4*A*D=1
%        N. Chernov and C. Lesort, "Least squares fitting of circles",
%        J. Math. Imag. Vision, Vol. 23, 239-251 (2005)
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)             
%             ParIni = [a b R] is the initial guess (supplied by user)
%     Output: Par = [a b R] is the fitting circle:                                               
%                           center (a,b) and radius R                                                                  
%                                                                                                     
    

                                                                                                                               
factorUp=10;  factorDown=0.04;     
lambda0=0.01;  epsilon=0.000001;    
IterMAX = 50;  AdjustMax = 20; 
Xshift=0;  Yshift=0;  dX=1;  dY=0;                                                                                    
n = size(XY,1);      # number of data points                                                                                   
                                                                                                                            
#     starting with the given initial guess
anew = ParIni(1) + Xshift;        
bnew = ParIni(2) + Yshift;                                                             
Anew = 1/(2*ParIni(3));                                                                                   
aabb = anew*anew + bnew*bnew;                                                                                  
Fnew = (aabb - ParIni(3)*ParIni(3))*Anew;                                                                
Tnew = acos(-anew/sqrt(aabb));                                                             
if (bnew > 0)                                                                                              
    Tnew = 2*pi - Tnew;                                              
end                                                                                                                  
%        if 1+4*Anew*Fnew < epsilon                                                                                     
%            fprintf(1,' +++ violation:  %f\n',1+4*Anew*Fnew);                                                                 
%        end                                                                                                                   
                                                
VarNew = VarCircle(XY,ParIni);                                                   
                                                                                                 
#     initializing lambda and iter                                                                   
                                                                                        
lambda = lambda0;  finish = 0;  
                                                                                                      
for iter=1:IterMAX                                                                               
                                                                      
    Aold = Anew;  Fold = Fnew;  Told = Tnew;  VarOld = VarNew;
                                                                                                                      
    H = sqrt(1+4*Aold*Fold);                                                                 
    aold = -H*cos(Told)/(Aold+Aold) - Xshift;
    bold = -H*sin(Told)/(Aold+Aold) - Yshift;
    Rold = 1/abs(Aold+Aold); 
#    fprintf(1,'%2d  (%f, %f)  %f  %.8f\n',iter,aold,bold,Rold,sqrt(VarOld));    
#           computing moments                                                                                         
         
    DD = 1 + 4*Aold*Fold; 
    D = sqrt(DD);  
    CT = cos(Told); 
    ST = sin(Told);    
    H11=0; H12=0; H13=0; H22=0; H23=0; H33=0; F1=0; F2=0; F3=0;           
                                                            
    for i=1:n                                                                                                
        Xi = XY(i,1) + Xshift;   
        Yi = XY(i,2) + Yshift;       
        Zi = Xi*Xi + Yi*Yi;  
        Ui = Xi*CT + Yi*ST;             
        Vi =-Xi*ST + Yi*CT;

        ADF = Aold*Zi + D*Ui + Fold;    
        SQ = sqrt(4*Aold*ADF + 1);           
        DEN = SQ + 1;                                          
        Gi = 2*ADF/DEN;   
        FACT = 2/DEN*(1 - Aold*Gi/SQ);      
        DGDAi = FACT*(Zi + 2*Fold*Ui/D) - Gi*Gi/SQ;                
        DGDFi = FACT*(2*Aold*Ui/D + 1);        
        DGDTi = FACT*D*Vi;    
                                                          
        H11 = H11 + DGDAi*DGDAi;                 
        H12 = H12 + DGDAi*DGDFi;                           
        H13 = H13 + DGDAi*DGDTi;                                          
        H22 = H22 + DGDFi*DGDFi;                  
        H23 = H23 + DGDFi*DGDTi;                                     
        H33 = H33 + DGDTi*DGDTi;                        
                                                 
        F1 = F1 + Gi*DGDAi;             
        F2 = F2 + Gi*DGDFi;    
        F3 = F3 + Gi*DGDTi;                                                        
    end                                                        

    for adjust=1:AdjustMax                                                               
                                                                              
%             Cholesly decomposition                                                                                           
                                                                                                                               
        G11 = sqrt(H11 + lambda);                                                                                       
        G12 = H12/G11;                                                                              
        G13 = H13/G11;                                                                                                 
        G22 = sqrt(H22 + lambda - G12*G12);                                                              
        G23 = (H23 - G12*G13)/G22;                                             
        G33 = sqrt(H33 + lambda - G13*G13 - G23*G23);                                                                          
                                                                                                                               
        D1 = F1/G11;                                                                                                           
        D2 = (F2 - G12*D1)/G22;                                                                                                
        D3 = (F3 - G13*D1 - G23*D2)/G33;                                                                                       

        dT = D3/G33;                                                                                                           
        dF = (D2 - G23*dT)/G22;                                                                                                
        dA = (D1 - G12*dF - G13*dT)/G11;                                                                                       
                                                                                                                               
%            updating the parameters
                                                                                                                               
        Anew = Aold - dA;                                                                                                      
        Fnew = Fold - dF;                                                                                                      
        Tnew = Told - dT;                                                                                                      
%        fprintf(1,'%2d   %.8f\n',iter,lambda);                                                                                

        if (1+4*Anew*Fnew < epsilon && lambda>1)                                                                               
%            fprintf(1,'     violation:  %f\n',1+4*Anew*Fnew);                                                                 
            Xshift = Xshift + dX;                                                                                              
            Yshift = Yshift + dY;                                                                                              
                                                                                                                               
            H = sqrt(1+4*Aold*Fold);                                                                                           
            aTemp = -H*cos(Told)/(Aold+Aold) + dX;                                                                             
            bTemp = -H*sin(Told)/(Aold+Aold) + dY;                                                                             
            rTemp = 1/abs(Aold+Aold);                                                                                          
                                                                                                                               
            Anew = 1/(rTemp + rTemp);                                                                                          
            aabb = aTemp*aTemp + bTemp*bTemp;                                                                                  
            Fnew = (aabb - rTemp*rTemp)*Anew;                                                                                  
            Tnew = acos(-aTemp/sqrt(aabb));                                                                                    
            if bTemp > 0                                                                                                       
               Tnew = 2*pi - Tnew;                                                                                             
            end                                                                                                                
            VarNew = VarOld;                                                                                                   
            break;                                                                                                             
        end                                                                                                                    
                                                                                                                               
        if 1+4*Anew*Fnew < epsilon                                                                                             
           lambda = lambda * factorUp;                                                                                         
           continue;                                                                                                           
        end                                                                                                                    

       DD = 1 + 4*Anew*Fnew;                                                                                                  
        D = sqrt(DD);                                                                                                          
        CT = cos(Tnew);                                                                                                        
        ST = sin(Tnew);                                                                                                        
                                                                                                                               
        GG = 0;                                                                                                                
                                                                                                                               
        for i=1:n                                                                                                              
            Xi = XY(i,1) + Xshift;                                                                                             
            Yi = XY(i,2) + Yshift;                                                                                             
            Zi = Xi*Xi + Yi*Yi;                                                                                                
            Ui = Xi*CT + Yi*ST;                                                                                                
                                                                                                                               
            ADF = Anew*Zi + D*Ui + Fnew;                                                                                       
            SQ = sqrt(4*Anew*ADF + 1);                                                                                         
            DEN = SQ + 1;                                                                                                      
            Gi = 2*ADF/DEN;                                                                                                    
            GG = GG + Gi*Gi;                                                                                                   
        end                                                                                                                    
                                                                                                                               
        VarNew = GG/(n-3);                                                                                                     
                                                                                                                               
        H = sqrt(1+4*Anew*Fnew);                                                                                               
        anew = -H*cos(Tnew)/(Anew+Anew) - Xshift;                                                                              
        bnew = -H*sin(Tnew)/(Anew+Anew) - Yshift;                                                                              
        Rnew = 1/abs(Anew+Anew);                                                                                               
                                                                                                                               
%             checking if improvement is gained                                                                                

        if VarNew <= VarOld      %   yes, improvement                                                                          
           progress = (abs(anew-aold) + abs(bnew-bold) + abs(Rnew-Rold))/(Rnew+Rold);                                          
           if progress < epsilon                                                                                               
              Aold = Anew;                                                                                                     
              Fold = Fnew;                                                                                                     
              Told = Tnew;                                                                                                     
              VarOld = VarNew; %#ok<NASGU>                                                                                     
              finish = 1;                                                                                                      
              break;                                                                                                           
           end                                                                                                                 
           lambda = lambda * factorDown;                                                                                       
           break;                                                                                                              
        else                     %   no improvement                                                                            
           lambda = lambda * factorUp;                                                                                         
           continue;                                                                                                           
        end                                                                                                                    
    end                                                                                                                        
    if finish == 1                                                                                                             
       break;                                                                                                                  
    end                                                                                                                        
end                                                                                                                       

H = sqrt(1+4*Aold*Fold);                                                                                        
Par(1) = -H*cos(Told)/(Aold+Aold) - Xshift;                                                      
Par(2) = -H*sin(Told)/(Aold+Aold) - Yshift;                                                 
Par(3) = 1/abs(Aold+Aold);                                                                                                     
                                                                                                                               
end  % LMA                                                                                                                     
"""
