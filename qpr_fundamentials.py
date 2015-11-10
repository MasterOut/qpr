# -*- coding: utf-8 -*-


import numpy as np

def norm2(v):       # calculate the euclidic norm
    v=np.power(v, 2)
    vsum = v.sum()
    return vsum

def fl(v):          # calculate the floor and return as integer value
    return np.floor(v).astype(int)

def nf(K, v):       # norm2 of floor (K*v)
    return norm2(fl(K * v))
    