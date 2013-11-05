#! /usr/bin/env python

import numpy
from numpy import linalg
import math
#import copy
#import spd


def get_zdata_segment(zdata, t_Arai, tmin, tmax):  # 
#    print "tmin, tmax", tmin, tmax
    lower_bound = t_Arai.index(tmin)
    upper_bound = t_Arai.index(tmax)
    zdata_segment = zdata[lower_bound:upper_bound + 1] # inclusive of upper bound
#    print "zdata segment from get_zdata_segment"
#    print zdata_segment
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
    return {'X1_prime': X1_prime, 'X2_prime': X2_prime, 'X3_prime': X3_prime, 'X1': X1, 'X2': X2, 'X3': X3, 'X1_mean': means[0], 'X2_mean': means[1], 'X3_mean': means[2], 'zdata_mass_center': means }

def get_orientation_tensor(X1_p, X2_p, X3_p):  # TRY DOING IT RON STYLE, with cov
    X1_p, X2_p, X3_p = numpy.array(X1_p), numpy.array(X2_p), numpy.array(X3_p)
    orient_tensor = [[sum(X1_p * X1_p), sum(X1_p * X2_p), sum(X1_p * X3_p)],
                     [sum(X1_p * X2_p), sum(X2_p * X2_p), sum(X2_p * X3_p)],
                     [sum(X1_p * X3_p), sum(X2_p * X3_p), sum(X3_p * X3_p)]]
    tau, V = numpy.linalg.eig(orient_tensor) 
# ATTEND!!  POSSIBLY THIS SHOULD BE numpy.linalg.eig(numpy.cov(orient_tensor))  numpy.cov gets the covariance matrix of data
    return {'orient_tensor': orient_tensor, 'tau': tau, 'V': V}


def tauV(T):
    """      
    gets the eigenvalues (tau) and eigenvectors (V) from matrix T
    """
    t,V,tr=[],[],0.
    ind1,ind2,ind3=0,1,2
    # tau, V
    evalues,evectmps=numpy.linalg.eig(T)
    evectors=numpy.transpose(evectmps)  # to make compatible with Numeric convention
    for tau in evalues:
        tr+=tau  # tr totals tau values
    if tr!=0:
        for i in range(3):
#            print "original:", evalues[i]
            evalues[i]=evalues[i] /tr  # convention is norming eigenvalues so they sum to 1.  why?  just because.  :)
#            print "normed:", evalues[i]
    else:
        return t,V  # if eigenvalues add up to zero, no sorting is needed
# sort evalues,evectors                                            
    t1,t2,t3=0.,0.,1.
    for k in range(3):  
        if evalues[k] > t1:
            t1,ind1=evalues[k],k
        if evalues[k] < t3:
            t3,ind3=evalues[k],k
#    print "ind1, ind2, ind3", ind1, ind2, ind3
#    print "t1, t2, t3", t1, t2, t3
    for k in range(3):
        if evalues[k] != t1 and evalues[k] != t3:
            t2,ind2=evalues[k],k
#    print "ind1, ind2, ind3", ind1, ind2, ind3
#    print "t1, t2, t3", t1, t2, t3
    V.append(evectors[ind1])
    V.append(evectors[ind2])
    V.append(evectors[ind3])
    t.append(t1)
    t.append(t2)
    t.append(t3)
#    print "tau", t
#    print "V", V
    return t,V


def sort_eigenvectors_and_eigenvalues(orient_tensor):
    # may need to sort eigenvectors, eigenvalues.  possibly not in order...?
    pass


# not sure if this is needed.  there is a method for getting direction along V1 using the vector dot product of V1 and R.  I just don't know if we need this.  
#def get_reference_vector(X1_prime, X2_prime, X3_prime):
#    n = len(X1_prime) - 1
#    X1 = X1_prime[0] - X1_prime[n]
#    X2 = X2_prime[0] - X2_prime[n]
#    X3 = X3_prime[0] - X3_prime[n]
#    R= numpy.array([X1, X2, X3])
#    return R
#    return numpy.array([0, 0, 0])
    
#tau, V = tauV(T['orient_tensor'])
#V1 = V[0]
#R = get_reference_vector(t['X1_prime'], t['X2_prime'], t['X3_prime'])
#dot = numpy.dot(V1, R) # dot product of reference vector and the principal axis of the V matrix
#if dot < -1:
#    dot = -1
#elif dot > 1:
#    dot = 1



ex_zdata = [[1., 1., 1.],[1., 2., 3.], [2., 4., 4.5], [2., 3., 7.5]]
real_zdata = [[-0.0178569 , -0.05399429,  0.9863136 ],
       [-0.02552379, -0.05074841,  0.95614538],
       [-0.02060363, -0.03101095,  0.96916464],
       [-0.05547474, -0.0306237 ,  0.979881  ],
       [-0.02034083, -0.03715329,  0.93277372],
       [-0.03429246, -0.03023287,  0.90245725],
       [-0.02914187, -0.03293889,  0.81206295]]


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
    Rs=numpy.sqrt(Xs**2+Ys**2+Zs**2) # calculate resultant vector length                 
    Decs=(numpy.arctan2(Ys,Xs)/rad)%360. # calculate declination taking care of correct quadrants (arctan2) and making modulo 360.
    try:
        Incs=numpy.arcsin(Zs/Rs)/rad # calculate inclination (converting to degrees) #                             
    except:
        print 'trouble in cart2dir' # most likely division by zero somewhere
        return numpy.zeros(3)
    return numpy.array([Decs,Incs,Rs]).transpose() # return the directions list

# Lisa sequence:
#T=numpy.array(Tmatrix(X))


def get_dec_and_inc(zdata, t_Arai, tmin, tmax, anchored=True):
    zdata_segment = get_zdata_segment(zdata, t_Arai, tmin, tmax)
    data = get_cart_primes_and_means(zdata_segment, anchored)
    means = data['zdata_mass_center']
    T = get_orientation_tensor(data['X1_prime'], data['X2_prime'], data['X3_prime'])
    tau,V=tauV(T['orient_tensor'])
    PDir=cart2dir(V[0])
    best_fit_vector = V[0]
    vector = adjust_best_fit_vector(best_fit_vector, means, zdata) # means is mass 
    if PDir[1] < 0:  # this whole transformatio is NOT done in thellier_gui.  ask Ron
        PDir[0]+=180. 
        PDir[1]=-PDir[1]
    PDir[0]=PDir[0]%360.
    dec = PDir[0]
    inc = PDir[1]
#    print "tau", tau
#    print "V", V
#    print "anchored:", anchored, " best fit vector", best_fit_vector, "adjusted vector", vector
    return dec, inc, vector, tau, V, means

def adjust_best_fit_vector(vector, cm, zdata):
    """finds correct sign "+/-" for best fit vector"""
#cm = zdata_mass_center
    cm = numpy.array(cm)
    best_fit_vector = numpy.array(vector)
    v1_plus = best_fit_vector * numpy.sqrt(sum(cm**2))
    v1_minus = best_fit_vector * -1. * numpy.sqrt(sum(cm**2))
    test_v = zdata[0] - zdata[-1]
#    print "test_v", test_v
    if numpy.sqrt(sum((v1_minus-test_v)**2)) < numpy.sqrt(sum((v1_plus-test_v)**2)):
        best_fit_vector = best_fit_vector* -1. 
    return best_fit_vector

#        cm=array(mean(zdata_segment.T,axis=1)) # center of mass  
 #       v1_plus=v1*sqrt(sum(cm**2))
#        v1_minus=v1*-1*sqrt(sum(cm**2))
#        test_v=zdata_segment[0]-zdata_segment[-1]
#        if sqrt(sum((v1_minus-test_v)**2)) < sqrt(sum((v1_plus-test_v)**2)):
#         DIR_PCA=self.cart2dir(v1*-1)
#         best_fit_vector=v1*-1
#        else:
#         DIR_PCA=self.cart2dir(v1)
#         best_fit_vector=v1



def get_MAD(tau): # works
    # tau is ordered so that tau[0] > tau[1] > tau[2]
# makes no difference whether use math.degrees or / rad to turn answer into degrees
#    rad = numpy.pi / 180.
    MAD = math.degrees(numpy.arctan(numpy.sqrt((tau[1] + tau[2]) / tau[0])) )# / rad # should work, AND DOES!!
    return MAD

def Lisa_get_MAD(t): # also works
#    print "Lisa_get_MAD tau:", t
    rad = numpy.pi / 180.
    s1=numpy.sqrt(t[0])
    MAD=numpy.arctan(numpy.sqrt(t[1]+t[2])/s1)/rad
#    print "Lisa MAD", MAD
    return MAD


# cm same as ...pars['zdata_mass_center']
def Lisa_get_DANG(cm, Dir): # working
   # print "cm, Dir", cm, Dir
    CMdir=cart2dir(cm)
    Dirp=cart2dir([Dir[0],Dir[1], Dir[2]])
    dang=pmag_angle(CMdir,Dirp) 
#    print "Lisa cmdir, dirp", CMdir, Dirp
    return dang

def Ron_get_DANG(cm, best_fit_vector):
    """gets deviation angle between center of mass and the free-float best-fit vector"""
    cmdir = cm
    dirp = best_fit_vector
    DANG=math.degrees( numpy.arccos( ( numpy.dot(cmdir, dirp) )/( numpy.sqrt(sum(cmdir**2)) * numpy.sqrt(sum(dirp**2)))))
#    print "Ron cmdir, dirp,", cmdir, dirp
    return DANG


def dir2cart(d): # from pmag.py
   # converts list or array of vector directions, in degrees, to array of cartesian coordinates, in x,y,z    
    ints=numpy.ones(len(d)).transpose() # get an array of ones to plug into dec,inc pairs             
    d=numpy.array(d)
    rad=numpy.pi/180.
    if len(d.shape)>1: # array of vectors                                                         
        decs,incs=d[:,0]*rad,d[:,1]*rad
        if d.shape[1]==3: ints=d[:,2] # take the given lengths                  
    else: # single vector                                 
        decs,incs=numpy.array(d[0])*rad,numpy.array(d[1])*rad
        if len(d)==3:
            ints=numpy.array(d[2])
        else:
            ints=numpy.array([1.])
    cart= numpy.array([ints*numpy.cos(decs)*numpy.cos(incs),ints*numpy.sin(decs)*numpy.cos(incs),ints*numpy.sin(incs)]).transpose()
    return cart


def pmag_angle(D1,D2): # use this 
    """          
    finds the angle between lists of two directions D1,D2
    """
    D1=numpy.array(D1)
    if len(D1.shape)>1:
        D1=D1[:,0:2] # strip off intensity
    else: D1=D1[:2]
    D2=numpy.array(D2)
    if len(D2.shape)>1:
        D2=D2[:,0:2] # strip off intensity
    else: D2=D2[:2]
#    print "D1", D1
#    print "D2", D2
    X1=dir2cart(D1) # convert to cartesian from polar
    X2=dir2cart(D2)
    angles=[] # set up a list for angles
    for k in range(X1.shape[0]): # single vector
        angle= numpy.arccos(numpy.dot(X1[k],X2[k]))*180./numpy.pi # take the dot product
#    angle=numpy.arccos( ( numpy.dot(v1, v2) )/( numpy.sqrt(math.fsum(v1**2)) * numpy.sqrt(math.fsum(v2**2))))    
        angle=angle%360.
        angles.append(angle)
    return numpy.array(angles)


def get_angle_difference(v1, v2):
    """returns angular difference in degrees between two vectors.  takes in cartesian coordinates."""
    v1 = numpy.array(v1)
    v2 = numpy.array(v2)
    print "v1:", v1, "v2:", v2
    print "numerator", numpy.dot(v1, v2)
    print "v1**2", v1**2
    print "sum(v1**2)", math.fsum(v1**2)
    print "numpy.sqrt(sum(v1**2))", numpy.sqrt(math.fsum(v1**2))
    print " * numpy.sqrt(sum(v2**2))", numpy.sqrt(math.fsum(v2**2))
    print "denominator", ( numpy.sqrt(math.fsum(v1**2)) * numpy.sqrt(math.fsum(v2**2)))
    angle=numpy.arccos( ( numpy.dot(v1, v2) )/( numpy.sqrt(math.fsum(v1**2)) * numpy.sqrt(math.fsum(v2**2))))    
    print "angle in radians:", angle
    return math.degrees(angle)


def get_alpha(anc_fit, free_fit): # 
    """returns angle between anchored best fit vector and free best fit vector"""
    alpha_deg = get_angle_difference(anc_fit, free_fit)
    return alpha_deg

def get_DANG(mass_center, free_best_fit_vector): # not working in all cases!  but working in some
    DANG = get_angle_difference(mass_center, free_best_fit_vector)
    return DANG


def get_NRM_dev(dang, x_avg, y_int):
    NRM_dev = ( numpy.sin(dang) * numpy.sqrt(x_avg[0]**2 + x_avg[1]**2 + x_avg[2]**2) ) / y_int
    NRM_dev *= 100.
    return NRM_dev

def get_theta(B_lab_dir, ChRM): # FINISH MEEEEEE
    B_lab_cart = dir2cart(B_lab_dir)
    get_angle_difference(B_lab_cart, cart) # you should change it so that get_angle_difference can take dir or cart
    pass


dir1 = [0, 90, 1]
cart1 = dir2cart(dir1)

dir2 = [1,2,3]
cart2 = dir2cart(dir2)


def get_gamma(B_lab_dir, pTRM_dir):
    B_lab_cart = dir2cart(B_lab_dir)
    pTRM_cart = dir2cart(pTRM_dir)
    print "B_lab_dir", B_lab_dir, "pTRM_dir", pTRM_dir
    print "B_lab_cart", B_lab_cart # [ ]
    print "pTRM cart", pTRM_cart # [[ ]] 
    gamma1 = pmag_angle(B_lab_dir, pTRM_dir)
    gamma2 = get_angle_difference(B_lab_cart, pTRM_cart) # problem is likely because of the zero value in B_lab
    # pmag_angle and get_angle difference are equivalent with non zero values
    gamma2 = numpy.array([gamma2])
    boo = gamma1 == gamma2
    print "boo", boo   # why the fuck is this false?  what the fuck is the difference between these things?
    print type(boo)
    if str(gamma1) == str(gamma2):
        print "success"
        return gamma1, gamma2
    else:
        print type(gamma1), type(gamma2)
        print type(gamma1[0]), type(gamma2[0])
        print gamma1.ndim, gamma2.ndim
        print gamma1.shape, gamma2.shape
        print gamma1, gamma2
        return gamma1  # or return gamma 2.  both work with unittests, so why the fuck don't they come out as equal, hmm?




#means = list(numpy.mean(zdata.T,axis=1))
#m=array(mean(CART_pTRMS_orig.T,axis=1)) 
#        v1_plus=v1*sqrt(sum(cm**2))
#        v1_minus=v1*-1*sqrt(sum(cm**2))
#        test_v=CART_pTRMS_orig[0]-CART_pTRMS_orig[-1]
#        if sqrt(sum((v1_minus-test_v)**2)) > sqrt(sum((v1_plus-test_v)**2)):
#         DIR_PCA=self.cart2dir(v1*-1)
#         best_fit_vector=v1*-1
#        else:
#         DIR_PCA=self.cart2dir(v1)
#         best_fit_vector=v1


# put this into get dec_and_inc
#cm = zdata_mass_center
#v1_plus = best_fit_vector * numpy.sqrt(sum(cm**2))
#v1_minus = best_fit_vector * -1. * numpy.sqrt(sum(cm**2))
#test_v = zdata[0] - zdata[-1]
#        if sqrt(sum((v1_minus-test_v)**2)) > sqrt(sum((v1_plus-test_v)**2)):
#          best_fit_vector = v1* -1.
#        else:
#          best_fit_vector = v1


if False:
    import spd
    thing = spd.thing
#    thing1 = spd.thing1
#    thing1.calculate_all_statistics()
    cm = thing.pars['zdata_mass_center']
    Dir = thing.pars['best_fit_vector_Free']
    r = Lisa_get_DANG(cm, Dir)
#    r2 = Ron_get_DANG(numpy.array(cm), numpy.array(Dir))
    r2 = Ron_get_DANG(numpy.array(cm), Dir)
    r3 = get_DANG(cm, Dir)
    
    # best fit vector is coming out wrong.  positive when should be negative
    print "Dir", Dir
    print "cm", cm
    print "which thing??", thing.s
    print "best_fit_vector", thing.pars['best_fit_vector_Free']
    print "Lisa r,", r
    print "Ron r", r2
    print "Lori r", r3
#    r1 = Ron_get_DANG(cm, Dir)
#    print "Ron DANG", r1
#def Lisa_get_DANG(cm, Dir):
