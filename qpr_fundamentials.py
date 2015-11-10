# -*- coding: utf-8 -*-


import numpy as np

def norm2(v):
    """
    Calculates the euclidic norm 2 of a vector v to the power of 2.
    """
    v=np.power(v, 2)
    vsum = v.sum()
    return vsum

def fl(v):
    """
    Calculates the floor and returnes as integer value
    """
    return np.floor(v).astype(int)

def nf(K, v):
    """
    Calculates the norm2 to the power of 2 of floor (K*v)
    """
    return norm2(fl(K * v))
    