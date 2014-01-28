#! /usr/bin/env python

import numpy
import math
#import spd

def get_zdata_segment(zdata, t_Arai, tmin, tmax):  # 
    lower_bound = t_Arai.index(tmin)
    upper_bound = t_Arai.index(tmax)
    zdata_segment = zdata[lower_bound:upper_bound + 1] # inclusive of upper bound
    return zdata_segment

#       get center of mass for principal components (

#def get_organized(X):# from pmag.py
#    cm = [0., 0., 0.]
#    print "X", X
#    for cart in X:
#        print "cm", cm
#        print "cart", cart
#        for l in range(3):
#            cm[l]+=cart[l]/len(X) # this division may be wrong
#    for k in range(len(X)):
#        for l in range(3):
#            X[k][l]=X[k][l]-cm[l]
#    print "final cm", cm
#    print "final X", X
#    return list(X), cm

#def Tmatrix(X): 
#    """      
#    gets the orientation matrix (T) from data in X      
#    """
#    T=[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
#    for row in X:
#        for k in range(3):
#            for l in range(3):
#                T[k][l] += row[k]*row[l]
#    return T

def get_cart_primes_and_means(zdata_segment, anchored=True):
# seek simpler syntax
    zdata = numpy.array(zdata_segment)
    X1 = zdata[:,0]
    X2 = zdata[:,1]
    X3 = zdata[:,2]
    means = list(numpy.mean(zdata.T,axis=1))
    if anchored == False:
        X1_prime = list(X1 - means[0])#X1_mean)
        X2_prime = list(X2 - means[1])
        X3_prime = list(X3 - means[2])
    else:
        X1_prime, X2_prime, X3_prime = list(X1), list(X2), list(X3)
    return {'X1_prime': X1_prime, 'X2_prime': X2_prime, 'X3_prime': X3_prime, 
            'X1': X1, 'X2': X2, 'X3': X3, 'X1_mean': means[0], 'X2_mean': means[1], 
            'X3_mean': means[2], 'zdata_mass_center': means }

def get_orientation_tensor(X1_p, X2_p, X3_p):  
    X1_p, X2_p, X3_p = numpy.array(X1_p), numpy.array(X2_p), numpy.array(X3_p)
    orient_tensor = [[sum(X1_p * X1_p), sum(X1_p * X2_p), sum(X1_p * X3_p)],
                     [sum(X1_p * X2_p), sum(X2_p * X2_p), sum(X2_p * X3_p)],
                     [sum(X1_p * X3_p), sum(X2_p * X3_p), sum(X3_p * X3_p)]]
    tau, V = numpy.linalg.eig(orient_tensor) 
    return {'orient_tensor': orient_tensor, 'tau': tau, 'V': V}

def tauV(T):
    """      
    gets the eigenvalues (tau) and eigenvectors (V) from matrix T
    """
    t,V,tr=[],[],0.
    ind1,ind2,ind3=0,1,2
    evalues,evectmps=numpy.linalg.eig(T)
    evectors=numpy.transpose(evectmps)  # to make compatible with Numeric convention
    for tau in evalues:
        tr += tau  # tr totals tau values
    if tr != 0:
        for i in range(3):
            evalues[i]=evalues[i] / tr  # convention is norming eigenvalues so they sum to 1.
    else:
        return t,V  # if eigenvalues add up to zero, no sorting is needed
# sort evalues,evectors      
    t1, t2, t3 = 0., 0., 1.
    for k in range(3):  
        if evalues[k] > t1:
            t1,ind1 = evalues[k],k
        if evalues[k] < t3:
            t3,ind3 = evalues[k],k
    for k in range(3):
        if evalues[k] != t1 and evalues[k] != t3:
            t2,ind2=evalues[k],k
    V.append(evectors[ind1])
    V.append(evectors[ind2])
    V.append(evectors[ind3])
    t.append(t1)
    t.append(t2)
    t.append(t3)
    return t,V

def get_PD_direction(X1_prime, X2_prime, X3_prime, PD):
    """takes arrays of X1_prime, X2_prime, X3_prime, and the PD.  
    checks that the PD vector direction is correct"""
    n = len(X1_prime) - 1
    X1 = X1_prime[0] - X1_prime[n]
    X2 = X2_prime[0] - X2_prime[n]
    X3 = X3_prime[0] - X3_prime[n]
    R= numpy.array([X1, X2, X3])
    dot = numpy.dot(PD, R) # dot product of reference vector and the principal axis of the V matrix
    if dot < -1:
        dot = -1
    elif dot > 1:
        dot = 1
    if numpy.arccos(dot) > numpy.pi / 2:
        PD = -1. * numpy.array(PD)
    return PD

def cart2dir(cart):
    """
    converts a direction to cartesian coordinates
    """
    cart=numpy.array(cart)
    rad=numpy.pi/180. # constant to convert degrees to radians
    if len(cart.shape)>1:
        Xs,Ys,Zs=cart[:,0],cart[:,1],cart[:,2]
    else: #single vector
        Xs,Ys,Zs=cart[0],cart[1],cart[2]
    Rs=numpy.sqrt(Xs**2 + Ys**2 + Zs**2) # calculate resultant vector length                 
    Decs=(numpy.arctan2(Ys,Xs)/rad)%360. # calculate declination taking care of correct quadrants (arctan2) and making modulo 360.
    try:
        Incs=numpy.arcsin(Zs/Rs)/rad # calculate inclination (converting to degrees) #                             
    except:
        print 'trouble in cart2dir' # most likely division by zero somewhere
        return numpy.zeros(3)
    return numpy.array([Decs,Incs,Rs]).transpose() # return the directions list

def get_dec_and_inc(zdata, t_Arai, tmin, tmax, anchored=True):
    zdata_segment = get_zdata_segment(zdata, t_Arai, tmin, tmax)
    data = get_cart_primes_and_means(zdata_segment, anchored)
    means = data['zdata_mass_center']
    T = get_orientation_tensor(data['X1_prime'], data['X2_prime'], data['X3_prime'])
    tau,V=tauV(T['orient_tensor'])
    PD = get_PD_direction(data['X1_prime'], data['X2_prime'], data['X3_prime'], V[0]) # makes PD + or -
    PDir=cart2dir(PD) # PDir is direction
    vector = PD # best fit vector / ChRM is cartesian
    dec = PDir[0]
    inc = PDir[1]
    return dec, inc, vector, tau, V, means

def get_MAD(tau):
    # tau is ordered so that tau[0] > tau[1] > tau[2]
    MAD = math.degrees(numpy.arctan(numpy.sqrt((tau[1] + tau[2]) / tau[0])) )# / rad # should work, AND DOES!!
    return MAD

def Lisa_get_MAD(t):
    rad = numpy.pi / 180.
    s1=numpy.sqrt(t[0])
    MAD=numpy.arctan(numpy.sqrt(t[1] + t[2])/s1) / rad
    return MAD

# cm same as ...pars['zdata_mass_center']
def Lisa_get_DANG(cm, Dir): # working
    CMdir=cart2dir(cm)
    Dirp=cart2dir([Dir[0],Dir[1], Dir[2]])
    dang=pmag_angle(CMdir,Dirp) 
    return dang

def Ron_get_DANG(cm, best_fit_vector):
    """gets deviation angle between center of mass and the free-float best-fit vector"""
    cmdir = cm
    dirp = best_fit_vector
    DANG = math.degrees(numpy.arccos((numpy.dot(cmdir, dirp)) / 
            (numpy.sqrt(sum(cmdir**2)) * numpy.sqrt(sum(dirp**2)))))
    return DANG

def dir2cart(d): # from pmag.py
    """converts list or array of vector directions, in degrees, to array of cartesian coordinates, in x,y,z form """
    ints = numpy.ones(len(d)).transpose() # get an array of ones to plug into dec,inc pairs             
    d = numpy.array(d)
    rad = numpy.pi / 180.
    if len(d.shape) > 1: # array of vectors                                                         
        decs, incs = d[:,0] * rad, d[:,1] * rad
        if d.shape[1] == 3: ints = d[:,2] # take the given lengths                  
    else: # single vector                                 
        decs, incs = numpy.array(d[0]) * rad, numpy.array(d[1]) * rad
        if len(d) == 3:
            ints = numpy.array(d[2])
        else:
            ints = numpy.array([1.])
    cart = numpy.array([ints * numpy.cos(decs) * numpy.cos(incs),
                       ints * numpy.sin(decs) * numpy.cos(incs),
                       ints * numpy.sin(incs)
                     ]).transpose()
    return cart

def pmag_angle(D1,D2): # use this 
    """          
    finds the angle between lists of two directions D1,D2
    """
    D1 = numpy.array(D1)
    if len(D1.shape) > 1:
        D1 = D1[:,0:2] # strip off intensity
    else: D1 = D1[:2]
    D2 = numpy.array(D2)
    if len(D2.shape) > 1:
        D2 = D2[:,0:2] # strip off intensity
    else: D2 = D2[:2]
    X1 = dir2cart(D1) # convert to cartesian from polar
    X2 = dir2cart(D2)
    angles = [] # set up a list for angles
    for k in range(X1.shape[0]): # single vector
        angle = numpy.arccos(numpy.dot(X1[k],X2[k]))*180./numpy.pi # take the dot product
        angle = angle%360.
        angles.append(angle)
    return numpy.array(angles)

def new_get_angle_diff(v1,v2):
    v1 = numpy.array(v1)
    v2 = numpy.array(v2)
    angle = numpy.arctan2(numpy.linalg.norm(numpy.cross(v1, v2)), numpy.dot(v1, v2))
    return math.degrees(angle)
    

def get_angle_difference(v1, v2):
    """returns angular difference in degrees between two vectors.  takes in cartesian coordinates."""
    v1 = numpy.array(v1)
    v2 = numpy.array(v2)
    angle=numpy.arccos((numpy.dot(v1, v2) ) / (numpy.sqrt(math.fsum(v1**2)) * numpy.sqrt(math.fsum(v2**2))))    
    return math.degrees(angle)

def get_alpha(anc_fit, free_fit): # 
    """returns angle between anchored best fit vector and free best fit vector"""
    alpha_deg = get_angle_difference(anc_fit, free_fit)
    return alpha_deg

def get_DANG(mass_center, free_best_fit_vector): 
    DANG = get_angle_difference(mass_center, free_best_fit_vector)
    return DANG

def get_NRM_dev(dang, x_avg, y_int):
    NRM_dev = (numpy.sin(numpy.deg2rad(dang)) * numpy.linalg.norm(x_avg)) / y_int
    NRM_dev *= 100.
    return NRM_dev

def get_theta(B_lab_dir, ChRM_cart): 
    B_lab_cart = dir2cart(B_lab_dir)
    #    print "B_lab_cart", B_lab_cart
    #    print "ChRM", ChRM
    #print "b lab cart"
    #print B_lab_cart
#    B_lab_cart = [0., 0., -1.] # adding this from Greig's code
    #print "new b lab cart"
    #print B_lab_cart
    #print "ChRM"
    #print ChRM
    ChRM_dir = cart2dir(ChRM_cart)
   # print 'B_lab_dir', B_lab_dir
   # print 'ChRM (dir)', ChRM_dir
   # print "B_lab_cart:", B_lab_cart
   # print 'ChRM (cart):', ChRM_cart
    theta1 = get_angle_difference(B_lab_cart, ChRM_cart) # you should change it so that get_angle_difference can take dir or cart
    theta2 = pmag_angle(B_lab_dir, ChRM_dir)
    print 'theta from cart:',  theta1
    print 'theta from dir:', theta2
    return theta1

def get_gamma(B_lab_dir, pTRM_dir):
    B_lab_cart = dir2cart(B_lab_dir)
    pTRM_cart = dir2cart(pTRM_dir)
    print 'B_lab_dir', B_lab_dir
    print 'pTRM_dir', pTRM_dir
    print 'B_lab_cart', B_lab_cart
    print 'pTRM_cart', pTRM_cart
    gamma1 = pmag_angle(B_lab_dir, pTRM_dir)
    gamma2 = get_angle_difference(B_lab_cart, pTRM_cart) # problem is likely because of the zero value in B_lab
    # pmag_angle and get_angle difference are equivalent with non zero values
    gamma2 = numpy.array([gamma2])
    gamma3 = new_get_angle_diff(B_lab_cart, pTRM_cart)
    print "gamma1 (from dir):", gamma1
    print "gqmma2 (from cart):", gamma2
    print "gamma3 ( with atan from cart): ", gamma3
#    print "gamma3", gamma3
    if gamma1 - gamma2 <= .0000001: # checks that the two methods of getting gamma return approximately equal results
        return float(gamma1)
    else:
        print "not equal with different methods"
        print gamma1
        print gamma2
        return float(gamma2)
