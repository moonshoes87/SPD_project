#! /usr/bin/env python

import numpy
from numpy import linalg
#import copy
import spd


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
    if anchored:
        X1_prime = list(X1 - means[0])#X1_mean)
        X2_prime = list(X2 - means[1])
        X3_prime = list(X3 - means[2])
    else:
        X1_prime, X2_prime, X3_prime = list(X1), list(X2), list(X3)
    return {'X1_prime': X1_prime, 'X2_prime': X2_prime, 'X3_prime': X3_prime, 'X1': X1, 'X2': X2, 'X3': X3, 'X1_mean': means[0], 'X2_mean': means[1], 'X3_mean': means[2], 'means': means }

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
    print "tau", t
    print "V", V
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
    T = get_orientation_tensor(data['X1_prime'], data['X2_prime'], data['X3_prime'])
    t,V=tauV(T['orient_tensor'])
    PDir=cart2dir(V[0])
#    print "PDir", PDir 
    if PDir[1] < 0:  # this whole transformatio is NOT done in thellier_gui.  ask Ron
        PDir[0]+=180.
        PDir[1]=-PDir[1]
    PDir[0]=PDir[0]%360.
    dec = PDir[0]
    inc = PDir[1]
    print "tau", t
    print "V", V
    return dec, inc, t

def get_MAD(tau):
    MAD = numpy.arctan((numpy.sqrt(tau[1]**2 + tau[2]**2)) / tau[0])
    return MAD

def Lisa_get_MAD(t):
    s1=numpy.sqrt(t[0])
    MAD=numpy.arctan(numpy.sqrt(t[1]+t[2])/s1)/rad
    print "Lisa MAD", MAD

def Ron_get_MAD(tau):
    tau1, tau2, tau3 = tau[0], tau[1], tau[2]
    MAD=math.degrees(numpy.arctan(numpy.sqrt((t2+t3)/t1)))
    print "Ron MAD", MAD


#dec1, inc1 = get_dec_and_inc(spd.thing.zdata, spd.thing.t_Arai, spd.thing.tmin, spd.thing.tmax)
#dec2, inc2 = get_dec_and_inc(spd.thing1.zdata, spd.thing1.t_Arai, spd.thing1.tmin, spd.thing1.tmax)
#dec3, inc3 = get_dec_and_inc(spd.thing1.zdata, spd.thing1.t_Arai, spd.thing1.tmin, spd.thing1.tmax, anchored=False)
