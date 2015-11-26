# -*- coding: utf-8 -*-


import numpy as np


FLOAT_TYPES = [type(int()), type(np.float32()), type(np.float64())]
INT_TYPES = [type(int()), type(np.int32()), type(np.int64())]
"""
General mathematic functions
"""
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
    Calculates the Computation Rate as known as the function:
    R(h, a) = 1/2 * log2(1/ (norm2(a) - (P*(h*a)^2))/(1+ P*norm2^2(h)) )
    """

    ha_dotprod_pow2 = np.power(np.dot(h, a), 2)
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

"""
Write to file
"""
#class file_dump:
#    def __init__(self, filename, cmdl_print=True):
#        self.filename = filename    # store filename as a class variable
#        self.cmdl_print = cmdl_print  # if true, the cmd line will be printed too
#        self.f = open(self.filename, "w")   # open class file
#
#    def def_header(self, *args):
#        self.head_str = ""
#        self.line_str = ""
#        self.line_str_a = ""
#        self.print_config = args
##        prec_string = ":.2"        
#        
#        last_dist = 0
#        for tup in self.print_config:
#            dist = tup[1]
#            typ = tup[2]    # float, int, array, ... for precision
#            tabs = abs(dist - last_dist)    # nr of tabs to next column
#            tabs_str = tabs*"\t"
#            column_key = tup[0]
#    
#            if typ == 'float':
#                precision = ":.3"
#            else:
#                precision = ""
#           
#            if column_key == "a":
#                additional = " : {anr}"
#                self.line_str_a += "\t{a}" + additional
#            else:
#                additional = "" 
#            self.line_str += tabs_str + "{" +  column_key + precision + "}" + additional
#            self.head_str += tabs_str + column_key
#            self.line_str_a += tabs_str
#            last_dist = dist
#        
#        self.print_line(self.head_str)
##        print self.head_str
#        #print self.line_str
#        
#    def print_line(self, linestr):
#        if self.cmdl_print:
#            print linestr
#        self.f.write(linestr+"\n")
#        
#    def print_table_line(self, **kwargs):
#        a_len = len(kwargs["a"]) # number of different coefficient vectors
#        a_dict = kwargs.pop("a")    # deletes and returns item with keyword "a"
#        a_print_dict = {}
#        kwargs["a"] = a_dict[0][0]      # für die erste Zeile, a
#        kwargs["anr"] = a_dict[0][1]    # für die erste Zeile, Anzahl des Auftretens vom ersten a
#        string = ""
#        if a_len == 1: # only one coefficient vector -> only one line needs to be printed            
#            string = self.line_str.format(**kwargs)
#            string += "\n"
#            if self.cmdl_print:
#                print(string)
#        else:
#            string = self.line_str.format(**kwargs) + "\n"
#            for i in range(1, a_len):               
#                a_print_dict[i] = {"a": a_dict[i][0], "anr": a_dict[i][1]}
#                string += self.line_str_a.format(**a_print_dict[i]) +"\n"
#        string += "\n"        
#        self.f.write(string)
#    
#    def print_dump(self, ITER_RANGE, **kwargs):
#        current_line_dict = {}   
#        self.print_line(self.head_str)
#        for i in ITER_RANGE:
#            for key in kwargs.keys():
#                current_data = kwargs[key][i]
#                current_line_dict[key] = kwargs[key][i]
#            self.print_table_line(current_line_dict)
#        
#    def close(self):
#        self.f.close()
#    
#
#class file_dumper:
#    def __init__(self, filename, cmdl_print=True):
#        self.filename = filename
#        self.cmdl_print = cmdl_print
#        self.f = open(filename, "w")
#        
#    def def_columns(self, column_tuple_list):
#        """
#        column_tuple_list is a list of tuples in the form: (column name, position)
#        """
#        self.col_names = []
#        self.col_pos = []
#        
#        for col_tup in column_tuple_list:
#            name = col_tup[0]
#            pos = col_tup[1]
#            self.col_names.append(name)
#            self.col_pos.append(pos)
#            if name in self.col_names:
#                print "ERROR: Column name is already defined."
#            if pos in self.col_pos:
#                print "ERROR: Column position has already been set."
#        print self.col_names
#        print self.col_pos







