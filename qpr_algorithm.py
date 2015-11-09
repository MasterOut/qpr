# -*- coding: utf-8 -*-


import numpy as np
import qpr_fundamentials as fun

def calc_r(u, L):   # calculate r
    if (L==2):
        r = u[L-1] / (1-fun.norm2(u[0])) * u[0]
    else:
        r = u[L-1] / (1-fun.norm2(u[0:L-1])) * u[0:L-1]
    return r

def calc_aKreuz(L, r, aKreuz):
    if (L==2):
        aKreuz[0] = r
        aKreuz[1] = 1
    else:
        aKreuz[0:L-1] = r
        aKreuz[L-1] = 1
    return aKreuz

def qp_relax(h, P, Ku, L, d):
    
    aKreuz = np.zeros(L)
    K = 0   # K
    Kl = 0  # running K
    
    hAbs = np.absolute(h)
    t = np.ones(L)
    
    np.copysign(t, h, t)    # copies sign elementwise from h to t, output is t
    
    p = hAbs.argsort(axis=None, kind='quicksort', order=None) # indexes of sorted channel vector
    hSort = np.copy(hAbs[p])   # sorted channel vector    
    
    b = 1 + P * fun.norm2(h)
    
    u = np.sqrt(P/b) * hSort
    
    r = calc_r(u, L)
    aKreuz = calc_aKreuz(L, r, aKreuz)
    
    # DETERMINE K
    if fun.nf(Ku, aKreuz) < b:
        K = Ku
    else:
        Kl = 1
        while (Ku != (Kl+1)):
            K = fun.fl(0.5*(Ku+Kl))  
            if fun.nf(K, aKreuz) < b:
                Kl = K
            else:
                Ku = K

        K = Kl
    
    
    aSquare = np.zeros(L, dtype=np.int)
    aSquare[L-1] = 1
    
    aKreuz = qp_quantization(aKreuz, u, L, K)
    
    fmin = fun.norm2(aSquare) - np.power( (np.dot(aSquare, u)), 2)
    f = fun.norm2(aKreuz) - np.power(d, 2)
    
    if f<fmin:
        aSquare = np.copy(aKreuz)
        fmin = f
    
    for l in range(0, L-1, 1):      # recover coefficient vector
        aSquare[p[l]] = t[p[l]] * aSquare[l]
    
    compRate = 0.5 * np.log(1/fmin)   # computation rate
    return [compRate, fmin]

def qp_quantization(aKreuz, u, L, K):
    # Quantization of a real-valued vector aKreuz of length L
    if K == 1:   # prevent from having empty range at K=1
        k_range = [1]
    else: 
        k_range = range(1, K+1)
    
    for k in k_range:
        aKreuz = k * aKreuz
        d = np.dot(aKreuz, u)

        for l in range(0, L, 1):    # 1 -> L-1
            v =  aKreuz[l]
            aKreuz[l] = np.floor(aKreuz[l])
            d = d + (aKreuz[l] - v) * u[l]
            x = 2* aKreuz[l] - 2* d*u[l] + 1 - np.power(u[l], 2)

            if x<0:
                aKreuz[l] += 1
                d += u[l]
    return aKreuz