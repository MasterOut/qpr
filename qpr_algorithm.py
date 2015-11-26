# -*- coding: utf-8 -*-


import numpy as np
import qpr_fundamentials as fun

def calc_r(u, L):   # calculate r
    """
    Calculates r as:
    r = u[L] / (1 - norm2^2( u[L-1] ) ) * u[1:L-1]
    """
    if (L==2):  # 0, 1 -> L-1=1
        r = u[L-1] / (1-fun.norm2(u[0])) * u[0]
    else:       # 
        r = u[L-1] / (1-fun.norm2(u[0:L-1])) * u[0:L-1]
    return r

def init_aC(L, r):
    """
    Returns the initial aC vector filled with r, except the last entry filled with 1.
    """
    aC = np.zeros(L)
    if (L==2):
        aC[0] = r
        aC[1] = 1
    else:
        aC[0:L-1] = r
        aC[L-1] = 1
    return aC

def normalize_vector(P, b, h_abs_sorted):
    """ Returns the normalized vector of h_abs """
    sqrt_P_b = np.sqrt(P/b)
    return sqrt_P_b * h_abs_sorted

def calc_fmin(aQ, u):
    aQ_norm2 = fun.norm2(aQ)
    dotprod = np.dot(aQ, u)
    return aQ_norm2 - np.power(dotprod, 2)

def qp_relax(h, P, Ku, L):
    """
    Implementation of the proposed algorithm by [ZM15].
    Returns a tuple with the calculated computation rate and the integer valued coefficient vector aQ.
    
    Parameters
    ----------
    h: array-like
        channel coefficient vector
    P: float
        Power
    Ku: int
        upper bound for K (maximal possible value in aQ)
    L: int
        length of channel vector/ aQ (=number of senders in the system)
    """
    K   = 0   # K
    Kl  = 0  # running K
    d   = 0.0
    
    h_abs   = np.absolute(h)# calculates the absolute of h
    t       = np.ones(L)    # signs of original channel vector
    a       = np.ones(L)    # coefficient vector
    t       = np.copysign(t, h)# copies signs of h elementwise to t (+/-1)
    
    p = np.argsort(h_abs, kind='quicksort') # indexes of sorted channel vector
    h_abs_sorted = np.copy(h_abs[p])   # ascending ordered vector  
    
    b = 1 + P * fun.norm2(h)    # part of denominator in Computation Rate formula in [1]
    
    u = normalize_vector(P, b, h_abs_sorted) # normalized channel vector
    
    r = calc_r(u, L)    # 
    aC1 = init_aC(L, r) # normalized initial a
    aC = np.zeros(L)    # used to test k*aC1
    
    """ DETERMINE K """
    if fun.nf(Ku, aC1) < b:
        K = Ku
    else:
        Kl = 1
        while (Ku != (Kl+1)):
            K = fun.fl(0.5*(Ku+Kl))  
            if fun.nf(K, aC1) < b:
                Kl = K
            else:
                Ku = K

        K = Kl
    
    """ Quantization """
    aQ = np.zeros(L, dtype=np.int)
    aQ[L-1] = 1
    
    """ norm2(aQ) = 1 , (aQ)^T * u = u[L-1] """
    # initially the maximum comp rate achievable is at fmin
    # fmin is reached if only one element in coeff vector is 1
    fmin = 1 - np.power(u[L-1], 2)

    if K == 1:   # prevent from having empty range at K=1
        k_range = [1]
    else: 
        k_range = range(1, K+1)
    
    for k in k_range:
        aC = k * aC1
        d = np.dot(aC, u)
        
        """ QUESTION:   in paper: for l <-1 to L-1
                        is in python: range(0, L) ? or range(0, L-1) 
        alg     python
        1       0
        2       1
        3       2
        …       …
        L-1     L-2
        L       L-1     range(0, L)   = 0, 1, 2 … L-1       L=4:    0,1,2,3
                        range(0, L-1) = 0, 1, 2 … L-1-1     L=4:    0,1,2
                        """
        for l in range(0, L-1):    # 1 -> L-1
            v =  aC[l]
            aC[l] = np.floor( aC[l] )
            d += (aC[l] - v) * u[l]
            # calculate condition term
            x = (2*aC[l]) - 2*d*u[l] + 1 - u[l]**2

            if x<0:
                aC[l] += 1
                d += u[l]
#            print "l={}, aC={}, x<0={}".format(l, aC, x<0)
        f = fun.norm2(aC) - np.power(d, 2)
        
#        print "\tR={0}, Rmax={1}, Rref={2}".format(comp_rate(f), comp_rate(fmin), fun.comp_rate(h_abs_sorted, aC, P))
#        print "\tk={}, fmin={}, f={}, aQ={}".format(k, fmin, f, aQ)
        if f<fmin:  # the less fmin, the bigger compRate
            aQ = np.copy(aC)
            fmin = f
        
    # recover coefficient vector
    for l in range(0,L):
        sign = t[p[l]]
        a[p[l]] = sign * aQ[l]

    compRate = comp_rate(fmin)   # computation rate
#    print "R={1:.4}, a={0}".format(a, compRate)
    return (compRate, a)

def comp_rate(f):
    return 0.5 * np.log2(1/f)

if __name__ == "__main__":
    h = np.array([1.2,0.3,0.8,2.1])
    L = 4
    Ku = 4
    qp_relax(h, 100, Ku, L)


