#!/usr/bin/env python 

import numpy
#import spd
#import lib_arai_plot_statistics as lib_arai

#thing = spd.thing
#  def get_xy_array(x_segment, y_segment):
#xy_real = lib_arai.get_xy_array(thing.x_Arai_segment, thing.y_Arai_segment)
#xy_real = [[0.12182379795523131, 0.98795180722891562], [0.1876743181924827, 0.95783132530120485], [0.21711848591334654, 0.96987951807228912], [0.32091411888460553, 0.98192771084337338], [0.48291503240724426, 0.9337349397590361], [0.72423703304017273, 0.90361445783132521], [1.0313987625711432, 0.81325301204819278]]    

xy_real = [[0.12182379795523131, 0.98795180722891562], [0.1876743181924827, 0.95783132530120485], [0.21711848591334654, 0.96987951807228912], [0.32091411888460553, 0.98192771084337338], [0.48291503240724426, 0.9337349397590361], [0.72423703304017273, 0.90361445783132521], [1.0313987625711432, 0.81325301204819278]]    

x_real = numpy.array([ 0.1218238 ,  0.18767432,  0.21711849,  0.32091412,  0.48291503, 0.72423703,  1.03139876])
y_real = numpy.array([ 0.98795181,  0.95783133,  0.96987952,  0.98192771,  0.93373494, 0.90361446,  0.81325301])

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
#    return { 'a':a,'b': b, 'r': r } #, XY[:,0], XY[:,1]
    return a,b,r
    
#    return (-1 * numpy.transpose((A[1:2]))) / A[0] / 2 + centroid
    


xy = numpy.array([[1,3],[0,2], [3,7], [8,10], [9,12]])
par = [1,2,7]

def VarCircle(XY = xy, Par = par):  # must have at least 4 sets of xy points or else division by zero occurs
    """
    computing the sample variance of distances from data points (XY) to the circle Par = [a b R]
    """
#    print "varcircle"
#    print "XY", XY
    n = len(XY)
#    print "n", n
    if n < 4:
        raise Warning("Circle cannot be calculated with less than 4 data points.  Please include more data")
    Dx = XY[:,0] - Par[0]
    Dy = XY[:,1] - Par[1]
    D = numpy.sqrt(Dx * Dx + Dy * Dy) - Par[2]
    result = numpy.dot(D, D)/(n-3)
    print n, Dx, Dy
    print D
#    print result
#    print "done varcircle"
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
    %     Input:  XY[n,2] is the array of coordinates of n points x[i]=XY[i,0]), y[i]=XY[i,1]             
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


    VarLambda = lambda0;  
    finish = 0;  
                                                                                                      
#    for iter=1:IterMAX
    for it in range(0,IterMAX):
                                                                      
        Aold = Anew  
        Fold = Fnew
        Told = Tnew
        VarOld = VarNew

        H = numpy.sqrt(1+4*Aold*Fold);                                                                 
        aold = -H*numpy.cos(Told)/(Aold+Aold) - Xshift;
        bold = -H*numpy.sin(Told)/(Aold+Aold) - Yshift;
        Rold = 1/abs(Aold+Aold); 
#        print "H, aold, bold, Rold", H, aold, bold, Rold

        DD = 1 + 4*Aold*Fold; 
        D = numpy.sqrt(DD);  
        CT = numpy.cos(Told); 
        ST = numpy.sin(Told);    
        H11=0; 
        H12=0; 
        H13=0; 
        H22=0; 
        H23=0; 
        H33=0; 
        F1=0; 
        F2=0; 
        F3=0;           
                                                            
        for i in range(0,n):
            Xi = XY[i,0] + Xshift;   
            Yi = XY[i,1] + Yshift;       
            Zi = Xi*Xi + Yi*Yi;  
            Ui = Xi*CT + Yi*ST;             
            Vi =-Xi*ST + Yi*CT;

            ADF = Aold*Zi + D*Ui + Fold;    
            SQ = numpy.sqrt(4*Aold*ADF + 1);           
            DEN = SQ + 1;                                          
            Gi = 2*ADF/DEN;   
            FACT = 2/DEN*(1 - Aold*Gi/SQ);      
            DGDAi = FACT*(Zi + 2*Fold*Ui/D) - Gi*Gi/SQ;                
            DGDFi = FACT*(2*Aold*Ui/D + 1);    # checked 
            DGDTi = FACT*D*Vi;    
                                                          
            H11 = H11 + DGDAi*DGDAi;                 
            H12 = H12 + DGDAi*DGDFi;                           
            H13 = H13 + DGDAi*DGDTi;                                          
            H22 = H22 + DGDFi*DGDFi;     # checked             
            H23 = H23 + DGDFi*DGDTi;                                     
            H33 = H33 + DGDTi*DGDTi;                        
                                                 
            F1 = F1 + Gi*DGDAi; # checked
            F2 = F2 + Gi*DGDFi;    
            F3 = F3 + Gi*DGDTi;
#            print "DGDFi", DGDFi, "H22", H22, "F1", F1 #


        for adjust in range(1,AdjustMax):
                                              
#             Cholesly decomposition                                     
                                                                       
            G11 = numpy.sqrt(H11 + VarLambda);
            G12 = H12/G11                                                                              
            G13 = H13/G11
            G22 = numpy.sqrt(H22 + VarLambda - G12*G12);                                                              
            G23 = (H23 - G12*G13)/G22;                                             
            G33 = numpy.sqrt(H33 + VarLambda - G13*G13 - G23*G23);                
                                                                                   
            D1 = F1/G11;                                            
            D2 = (F2 - G12*D1)/G22;                                                              
            D3 = (F3 - G13*D1 - G23*D2)/G33;                

            dT = D3/G33;  # checked                                  
            dF = (D2 - G23*dT)/G22 # checked
            dA = (D1 - G12*dF - G13*dT)/G11 # checked    
#            print 'dT', dT, 'dF', dF, 'dA', dA 
                                                                                   
#            updating the parameters
                                                                                            
            Anew = Aold - dA;  # checked
            Fnew = Fold - dF;                             
            Tnew = Told - dT;
#            print "Anew", Anew


            if 1+4*Anew*Fnew < epsilon and VarLambda>1:  # not hitting this condition in example or my code
                print "1+4*Anew*Fnew < epsilon and VarLambda>1:"
#            fprintf(1,'     violation:  %f\n',1+4*Anew*Fnew);             
                Xshift = Xshift + dX;                                          
                Yshift = Yshift + dY;                                                                               
                                                                                     
                H = numpy.sqrt(1+4*Aold*Fold);                               
                aTemp = -H*numpy.cos(Told)/(Aold+Aold) + dX;                                     
                bTemp = -H*numpy.sin(Told)/(Aold+Aold) + dY;                                      
                rTemp = 1/abs(Aold+Aold);                                       
                                                                             
                Anew = 1/(rTemp + rTemp);                         
                aabb = aTemp*aTemp + bTemp*bTemp;                          
                Fnew = (aabb - rTemp*rTemp)*Anew;                             
                Tnew = numpy.arccos(-aTemp/sqrt(aabb));                                       
                if bTemp > 0:
                    Tnew = 2*pi - Tnew;           
                VarNew = VarOld;                                         
                break;                               
           # end                                                    
            
            if 1+4*Anew*Fnew < epsilon:  # not hitting this condition in example or my code
                print "it", it
                VarLambda = VarLambda * factorUp;             
                print "VarLambda", VarLambda
                continue;              
           # end

            DD = 1 + 4*Anew*Fnew;                  
            D = numpy.sqrt(DD);                                                         
            CT = numpy.cos(Tnew);                                
            ST = numpy.sin(Tnew);    
                    
            GG = 0;                
                            
#            for i=1:n                 
            for i in range(0, n):
                Xi = XY[i,0] + Xshift;          
                Yi = XY[i,1] + Yshift;    
                Zi = Xi*Xi + Yi*Yi; # checked
                Ui = Xi*CT + Yi*ST;            
                                                 
                ADF = Anew*Zi + D*Ui + Fnew;                    
                SQ = numpy.sqrt(4*Anew*ADF + 1);               
                DEN = SQ + 1;                   
                Gi = 2*ADF/DEN; # checked  
                GG = GG + Gi*Gi;
              #  print "Zi", Zi, "Gi", Gi
           # end              
                                   
            VarNew = GG/(n-3);    
         
            H = numpy.sqrt(1+4*Anew*Fnew);               
            anew = -H*numpy.cos(Tnew)/(Anew+Anew) - Xshift;  #checked
            bnew = -H*numpy.sin(Tnew)/(Anew+Anew) - Yshift;  #checked
            Rnew = 1/abs(Anew+Anew); 
            print "anew", anew, "bnew", bnew

            if VarNew <= VarOld: #   yes, improvement                
               # print "VarNew <= VarOld", VarNew, VarOld
                progress = (abs(anew-aold) + abs(bnew-bold) + abs(Rnew-Rold))/(Rnew+Rold);      
                if progress < epsilon: 
                    print "Progress < epsilon"
                    Aold = Anew;          
                    Fold = Fnew;      
                    Told = Tnew;           
                    VarOld = VarNew # %#ok<NASGU>  
                    finish = 1;     
                    break;  

                VarLambda = VarLambda * factorDown
                break  
            else:                 #    %   no improvement  
                print "doing else"
                VarLambda = VarLambda * factorUp;      
                continue;     

        if finish == 1:
            break

    H = numpy.sqrt(1+4*Aold*Fold);                                                                                        
    result_a = -H*numpy.cos(Told)/(Aold+Aold) - Xshift;                                                      
    result_b = -H*numpy.sin(Told)/(Aold+Aold) - Yshift;                                                 
    result_r = 1/abs(Aold+Aold);       

    print result_a, result_b, result_r
    return result_a, result_b, result_r


new_xy = numpy.array([[.4, 6], [.3, 5.5], [.51, 5], [.7, 4.2], [.3, 3], [.8, 2.1]])
new_par_ini = [58.4404, 8.2870, 58.0924]
#a, b, r = LMA(new_xy, new_par_ini)
#print "new results"
#print a, b, r


def AraiCurvature(x,y):
    """
% Function for calculating the radius of the best fit circle to a set of 
% x-y coordinates.
% Paterson, G. A., (2011), A simple test for the presence of multidomain
% behaviour during paleointensity experiments, J. Geophys. Res., in press,
% doi: 10.1029/2011JB008369
% input: x = [1, 2, 3], y = [6, 4, 3]
% output: curvature, a (x coordinate circle center), b (y coordinate circle center), SSE (goodness of fit)
% parameters[0] = k
% parameters[1] = a
% parameters[2] = b
% parameters[3] = SSE (goodness of fit)
    """

    XY = []
    for num, X in enumerate(x):
        XY.append([X, y[num]])
    XY = numpy.array(XY)
    # Normalizevectors
    XY[:,0] = XY[:,0] / max(XY[:,0]) # important that values are all floats  
    XY[:,1] = XY[:,1] / max(XY[:,1]) # norms y values  # ThIS MAY ALREADY BE DONE.... from the thellier_magic processing
    X = XY[:,0]
    Y = XY[:,1]
                  
    #Provide the intitial estimate
    E1=TaubinSVD(XY);

    print "E1", E1

    #Determine the iterative solution
    E2=LMA(XY, E1);

    estimates=[E2[2], E2[0], E2[1]];
    
    best_a = E2[0]
    best_b = E2[1]
    best_r = E2[2]

    print "mean x"
    print numpy.mean(X)

    if best_a <= numpy.mean(X) and best_b <= numpy.mean(Y):
        print "-1/r"
        k = -1./best_r
    else:
        print "1/r"
        k = 1./best_r

    SSE = get_SSE(best_a, best_b, best_r, X, Y)

    #end
    print "best_r", best_r
    print "k, best_a, best_b, SSE"
    print k, best_a, best_b, SSE
    return k, best_a, best_b, SSE


a, b, r = 1,2, 3
x = [0.1, 0.2, 0.5]
y = [6.0, 4.0, 3.0]

def get_SSE(a,b,r,x,y):
    SSE = 0
    X = numpy.array(x)
    Y = numpy.array(y)
    for i in range(len(X)):
        x = X[i]
        y = Y[i]
        v = (numpy.sqrt( (x -a)**2 + (y - b)**2 ) - r )**2
#            print v                                                                                  
        SSE += v
    print SSE
    return SSE



xy = [[1,1],[0.,-2.], [.5,3.], [4., 5.]]
x_arai = [1,0.,.5,4.]
y_arai = [1., -2., 4., 5.]
#AraiCurvature(x_arai, y_arai)

