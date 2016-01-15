# -*- coding: utf-8 -*-


import numpy as np


FLOAT_TYPES = [type(int()), type(np.float32()), type(np.float64())]
INT_TYPES = [type(int()), type(np.int32()), type(np.int64())]

def norm2(v):
    """ Calculates the euclidic norm 2 of a vector v to the power of 2. """
    v=np.power(v, 2)
    return v.sum()

def fl(v):
    """ Calculates the floor and returnes as integer value """
    return np.floor(v).astype(int)

def nf(K, v):
    """ Calculates the norm2 to the power of 2 of floor (K*v) """
    return norm2( fl(K * v) )

def log2_plus(v):
    """  Returns the log in respect to base 2, if this is bigger or equal 0, whereas the function returns 0 if it is smaller than 0. """
    log2_v = np.log2(v)
    if log2_v < 0:
        ret_val = 0
    else:
        ret_val = log2_v
    return ret_val

""" Algorithm specific functions """    
def choose_Ku(L):
    """ Returns a upper bound Ku as listed in Table 1 in [1] """
    res = 0
    L_Ku_dict = {2:2, 3:3, 4:4, 5:5, 6:5, 7:5, 8:6, 9:6, 10:6, 11:6, 12:7, 13:6, 14:6, 15:6, 16:4}
    if (L in range(2, 17, 1)):
        res = L_Ku_dict[L]
    else:
        res = 0
    return res

def qpr_lin_fit(x, m, n):
    """ Linear function model to fit measuremnt data.  """
    return x*m + n
        
def comp_rate(h, a, P):
    """
    Calculates the Computation Rate via function:
    R(h, a) = 1/2 * log2(1/ (norm2(a) - (P*(h*a)^2))/(1+ P*norm2^2(h)) )
    """
    h = np.absolute(h)
    a = np.absolute(a)
    ha_dot = np.dot(h, a)
    ha_dotprod_pow2 = np.power(ha_dot, 2)
    denominator = norm2(a) - ( P * ha_dotprod_pow2 / (1 + P * norm2(h) ) )
    R = 0.5 * log2_plus(1 / denominator)
    return R 

def comp_rate_reference(h, a):
    """
    Iterates over a logscaled input P and calculates the computation rate with given parameters.
    
    Parameters
    ----------
    h: np.array
        channel vector
    a: np.array
        integer valued coefficient vector
    """
    P_NUM = 50
    P = np.logspace(0, 2, num=P_NUM)
    PdB = 10* np.log10(P)
    compRate_arr = np.zeros(P_NUM)
    
    for p_idx, p_val in enumerate(P):
        compRate_arr[p_idx] = comp_rate(h, a, p_val)    
    return PdB, compRate_arr

def cnt_cofVec_occ(a, di="increasing"):
    """
    Returns a dict called p_dict which contains data in the following form:
    
    p_dict = {nr_entry: [ [a_vec], nr_occur ]}
    
    Also a sorted list with entry numbers is returned.
    
    Parameters
    ----------
    a: array like
        array of coefficient vectors, which should be counted
    di: string
        gives the direction of sorting, increasing or decreasing is possible
    """
    nr_entry = -1
    p_dict = {}    
    for i in range(0, len(a)):  # iteration over array a
        found = False
        ai = a[i]               # current value of ai 
        p_dict_keys = p_dict.keys()
        
        for k in p_dict_keys:
            p_dict_value = p_dict[k]    # is a list of [a_vec], nr_occurence
            array = p_dict_value[0]
            if np.array_equal(ai, array):
                p_dict_value[1] += 1    # increment number of occurence
                found = True
                break               
            
        if not found:   # if not found in current listed array -> new entry
            nr_entry += 1
            p_dict[nr_entry] = [a[i], 1]
        
    sort_vec = sort_cofVec_occ(p_dict, di)           
    return p_dict, sort_vec

def sort_cofVec_occ(p_dict, direction):
    """ Return a list with indices to sort coefficient vectors """        
    sort_vec = [None] * len(p_dict)   # will be sorted
    for k in p_dict.keys(): # get the values of occurence for each ćoeff vector
        sort_vec[k] = p_dict[k][1]   # elements in list have the same index, as key in p_dict !!!
        
    if direction=='decreasing':
        for i in range(0, len(sort_vec)):
            sort_vec[i] = -sort_vec[i]    # negate to have a decreasing order
            
    ind_vec = np.argsort(sort_vec, axis=0) # sort sort_vec and save the indexes
          
    return ind_vec
    
def cnt_appearance(a, di=True):
    """
    Returns two lists. 1: list with channel vectors. 2: list with occurence of each channel vector
        
    Parameters
    ----------
    a: array like
        array of coefficient vectors, which will be counted
    di: bool
        True: decreasing, False: increasing order
    """
    lis = [] # [ [a_vector, cnt(a_vector] ),  ...]
    lis_app = []
    for i in range(0, len(a)):  # for all coefficient vectors in a:
        found = False
        ai = a[i]
        #ai_neg = np.negative(ai)
        
        for i, e in enumerate(lis): # if current coeff vector already seen, increment occurence list
            equal = array_equal_pos_neg(e, ai)
            if equal in [1, -1]:
                lis_app[i] += 1
                found = True
                break
        
        if not found:   # if not seen, add coeff vector to list and 1 to occurence list
            lis.append(ai)
            lis_app.append(1)
    
    lis_app_unsort = np.copy(lis_app)
    if di:  # if decreasing order is wanted: negate occurences
        lis_app = np.negative(lis_app)
        
    lis_sort_ind = np.argsort(lis_app, axis=0) # get indizes to sort lis_app
    
    lis_unsort = np.copy(lis)
    
    for i in range(0, len(lis)):
        lis[i] = lis_unsort[lis_sort_ind[i]]    
        lis_app[i] = lis_app_unsort[lis_sort_ind[i]]
    
    return lis, lis_app


def array_equal_pos_neg(a, b):
    ret = 0
    if np.array_equal(a, b):
        ret = 1
    elif np.array_equal(a, np.negative(b)):
        ret = -1
    return ret

