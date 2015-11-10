# -*- coding: utf-8 -*-


import numpy as np
import qpr_fundamentials as fun

def calc_r(u, L):   # calculate r
    if (L==2):
        r = u[L-1] / (1-fun.norm2(u[0])) * u[0]
    else:
        r = u[L-1] / (1-fun.norm2(u[0:L-1])) * u[0:L-1]
    return r

def calc_aCross(L, r, aCross):
    if (L==2):
        aCross[0] = r
        aCross[1] = 1
    else:
        aCross[0:L-1] = r
        aCross[L-1] = 1
    return aCross

def qp_relax(h, P, Ku, L):
    
    aCross = np.zeros(L)
    K = 0   # K
    Kl = 0  # running K
    d = 0.0
    
    hAbs = np.absolute(h)
    t = np.ones(L)
    
    np.copysign(t, h, t)    # copies sign elementwise from h to t, output is t
    
    p = hAbs.argsort(axis=None, kind='quicksort', order=None) # indexes of sorted channel vector
    hSort = np.copy(hAbs[p])   # sorted channel vector    
    
    b = 1 + P * fun.norm2(h)
    
    u = np.sqrt(P/b) * hSort
    
    r = calc_r(u, L)
    aCross = calc_aCross(L, r, aCross)
    
    # DETERMINE K
    if fun.nf(Ku, aCross) < b:
        K = Ku
    else:
        Kl = 1
        while (Ku != (Kl+1)):
            K = fun.fl(0.5*(Ku+Kl))  
            if fun.nf(K, aCross) < b:
                Kl = K
            else:
                Ku = K

        K = Kl
    
    # Quantization
    aSquare = np.zeros(L, dtype=np.int)
    aSquare[L-1] = 1
    fmin = fun.norm2(aSquare) - np.power( (np.dot(aSquare, u) ), 2)
    
    if K == 1:   # prevent from having empty range at K=1
        k_range = [1]
    else: 
        k_range = range(1, K+1)
    
    for k in k_range:
        aCross = k * aCross
        d = np.dot(aCross, u)

        for l in range(0, L):    # 1 -> L-1
            v =  aCross[l]
            aCross[l] = np.floor( aCross[l] )
            d = d + (aCross[l] - v) * u[l]  # type(d) = float64
            x = 2* aCross[l] - 2* d*u[l] + 1 - np.power(u[l], 2)

            if x<0:
                aCross[l] += 1
                d += u[l]
                
        f = fun.norm2(aCross) - np.power(d, 2) 
        if f<fmin:
            aSquare = np.copy(aCross)
            fmin = f
  

    for l in range(0, L):      # recover coefficient vector
        aSquare[p[l]] = t[ p[l] ] * aSquare[l]

    compRate = 0.5 * np.log(1/fmin)   # computation rate

    return (compRate, aSquare)
