# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 20:21:45 2015

@author: johannes
"""
import numpy as np
import qpr_algorithm as alg
import qpr_fundamentials as fun
import time
import qpr_csv_dump as qcsv
#import matplotlib.pyplot as plt
#import csv


def qpr_main(par_dict):
    """
    cr: Computation Rate
    a: coefficient vector
    h: channel vector
    P: Power
    Ku: upper bound for L
    L: width of channel vector
    """
    date_str = time.strftime("%Y_%m_%d_%H%M%S")
#    date_str = ""
    filename = 'run_' + date_str + '.txt'
#    fieldnames = ['nr', 'P', 'Pdb', 'h', 'a', 'a_occ', 'R', 'time', 'R_ref']
    w = qcsv.csv_dict_writer(filename)    
    
    
    its = par_dict["its"]    
    
    P = np.logspace(par_dict["pstart"], par_dict["pend"], par_dict["pnum"])
    Pdb = 10* np.log10(P)
    
    if par_dict["h"] == np.array([None]):
        use_std_normal = True
    else:
        use_std_normal = True   # currently not implemented

    L_list = range(par_dict["Lstart"], par_dict["Lend"]+1, par_dict["Lstep"])

    time_scale = par_dict["time_scale"]    
    
    # dict with results to write in csv file  
    w_dict = {}
    
    # array to hold cr_ref results before averaging
    cr_ref  = np.zeros((its, 1))

    its_range = range(0, its)
    
    """ main loop running algorithm and measure time and plot processing results """
    for L in L_list:
        print "L={}".format(L)
        par_dict["L"] = L
        # choose upper bound
        Ku = fun.choose_Ku(L)
        par_dict["Ku"] = Ku # to print in csv file 
        # print settings to csv file
        w.write_parameters(par_dict)
        for ind, p in enumerate(P):   
    #        pdb = 10*np.log10(p)           
            if use_std_normal:
                h = np.random.standard_normal(size=(its, L))  
                
            cr = np.zeros((its,1))  # array to hold cr values
            a = np.zeros((its, L))  # array to hold a coefficients
    
            """ iteration loopt with time measurement"""
            tstart = time.time()
            for i in its_range:
                cr[i], a[i] =  alg.qp_relax(h[i], p, Ku, L)
            tend = time.time()
            time_needed = (tend - tstart) * time_scale
            
            a_list, a_occ_list = cnt_occurence(a)
            cr_av = np.average(cr)
                    
            """ Comparison with reference computation rate formula """
            a_most = a_list[0]
            for i_h, v_h in enumerate(h):
                cr_ref[i_h] = fun.comp_rate(np.absolute(v_h), a_most, p)
            cr_ref_av = np.average(cr_ref)
            
            """ fill write dictionary to dump in file"""        
            w_dict = {'nr': ind,'P': p,'Pdb': Pdb[ind],'h': ['std'],'a': a_list,'a_occ': a_occ_list,'R': cr_av,'time': time_needed,'R_ref': cr_ref_av,}
            
            """ file dump """
            w.write_row(w_dict)
        
    """ Preprocess """
#    r = qcsv.csv_reader(filename, fieldnames)
#    r.plot_r_p()


def cnt_occurence(a, di=True):
    """
    Returns a dict called p_dict which contains data in the following form:
    
    p_dict = {nr_entry: [ [a_vec], nr_occur ]}
    
    Also a sorted list with entry numbers is returned.
    
    Parameters
    ----------
    a: array like
        array of coefficient vectors, which should be counted
    di: bool
        True: decreasing, False: increasing order
    """
    lis = [] # [ [a_vector, cnt(a_vector] ),  ...]
    lis_occ = []
    for i in range(0, len(a)):  # for all coefficient vectors in a:
        found = False
        ai = a[i]
        
        for i, e in enumerate(lis): # if current coeff vector already seen, increment occurence list
            if np.array_equal(e, ai):
                lis_occ[i] += 1
                found = True
                break
        
        if not found:   # if not seen, add coeff vector to list and 1 to occurence list
            lis.append(ai)
            lis_occ.append(1)
    
    
    lis_occ_unsort = np.copy(lis_occ)
    if di:  # if decreasing order is wanted: negate occurences
        lis_occ = negative(lis_occ)
        
    lis_sort_ind = np.argsort(lis_occ, axis=0) # get indizes to sort lis_occ
    
    lis_unsort = np.copy(lis)
    
    for i in range(0, len(lis)):
        lis[i] = lis_unsort[lis_sort_ind[i]]    
        lis_occ[i] = lis_occ_unsort[lis_sort_ind[i]]
    
    return lis, lis_occ

#==============================================================================
#  MAIN
#==============================================================================
if __name__ == '__main__':
    print(40*"-"+"\nQuadratic Programming Relaxation Approach to Compute-and-Forward Network Coding Design\n"+40*"-")
        
    par_dict = {
        "pstart": 0,    # dB/10
        "pend": 2,      # dB/10
        "pnum": 200, # number of points between pstart and pend
#        "L": 4,
        "Lstart": 4,
        "Lend": 16,
        "Lstep": 2,
        "h": np.array([None]),  # when None is used -> standard normal distribution will be used to create h's
        "its": 1000, # number of iterations per setting (h will be computed newly if std norm is used)
        "time_scale": 100, # 1000 for ms, 1 for s ...
    }
    
    # run simulation
    qpr_main(par_dict)
    
"""
[Computation rate: Koeffizienten a]

h = [0.4,1.3,1.2,0.6] mit Leistung P = 100
---------------------
1.3044: (1, 2, 2, 1)    -> found
1.0447: (1, 4, 4, 2)
0.9729: (1, 3, 3, 1)
0.8871: (0, 1, 1, 0)
0.7248: (2, 6, 6, 3)
0.7088: (0, 1, 1, 1)

h = [1.2,0.3,0.8,2.1] mit Leistung P = 100
---------------------
1.3530: (3, 1, 2, 5)    -> found (with Ku=5 !)
1.2914: (1, 0, 1, 2)    -> found (with Ku=4 as said in table in paper)
1.1123: (4, 1, 3, 7)
0.9009: (3, 1, 2, 6)
0.8565: (2, 0, 1, 3)
0.8046: (5, 1, 3, 8)

h = [0.1,0.3,0.8,0.5] mit Leistung P = 100
---------------------
1.0294: (0, 1, 2, 1)    -> found
0.8448: (0, 1, 3, 2)
0.8448: (0, 0, 1, 1)    
0.7370: (0, 0, 1, 0)
0.5922: (0, 1, 1, 1)
0.4183: (1, 2, 5, 3)
""" 

















