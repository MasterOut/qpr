# -*- coding: utf-8 -*-


import numpy as np
import qpr_algorithm as alg
import qpr_fundamentials as fun
import time
import qpr_csv_dump as qcsv


def qpr_main(par_dict):
    """
    Runs qpr_relax in 2 loops over L (channel vector length) and P.
    par_dict holds settings for the simulation.
    """
    date_str = time.strftime("%Y_%m_%d_%H%M%S")
#    date_str = ""
    filename = 'run_' + date_str + '.txt'
#    fieldnames = ['nr', 'P', 'Pdb', 'h', 'a', 'a_occ', 'R', 'time', 'R_ref']
    w = qcsv.csv_dict_writer(filename)    
    
    # generate P and PdB
    P = np.logspace(par_dict["pstart"], par_dict["pend"], par_dict["pnum"])
    Pdb = 10* np.log10(P)
    
    # use standard normal distributed channel vector or given vector
    if np.size(par_dict["h"]) == 1:
        use_std_normal = True
        L_list = range(par_dict["Lstart"], par_dict["Lend"]+1, par_dict["Lstep"])
        its = par_dict["its"]
    else:
        use_std_normal = False   # currently not implemented
        L_list = [np.size(par_dict["h"])]
        its = 1
    
    # time is in seconds: multiply by time_scale to get ms (1000)
    time_scale = par_dict["time_scale"]    
    
    # dict with results to write in csv file  
    w_dict = {}
    
    # array to hold cr_ref results before averaging
    cr_ref  = np.zeros((its, 1))
    cr_ref_av = 0
    
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
            else:
                h = [ np.array(par_dict["h"]) ]
                
            cr = np.zeros((its,1))  # array to hold cr values
            a = np.zeros((its, L))  # array to hold a coefficients
    
            """ iteration loop with time measurement"""
            tstart = time.time()
            for i in its_range:
                cr[i], a[i] =  alg.qp_relax(h[i], p, Ku, L)
            tend = time.time()
            time_needed = (tend - tstart) * time_scale
            
            a_list, a_occ_list = fun.cnt_occurence(a)
            cr_av = np.average(cr)
            
            if par_dict["calc_ref"]:
                """ Comparison with reference computation rate formula """
                a_most = a_list[0]
                for i_h, v_h in enumerate(h):
                    cr_ref[i_h] = fun.comp_rate(np.absolute(v_h), a_most, p)
                cr_ref_av = np.average(cr_ref)
            
            # fill write_dict with results     
            w_dict = {'nr': ind,'P': p,'Pdb': Pdb[ind],'h': ['std'],'a': a_list,'a_occ': a_occ_list,'R': cr_av,'time': time_needed,'R_ref': cr_ref_av,}
            
            # write to file
            w.write_row(w_dict)
        """
        Preprocess with running qpr_csv_dump.py
        """

#==============================================================================
#  MAIN
#==============================================================================
if __name__ == '__main__':
    print(40*"-"+"\nQuadratic Programming Relaxation Approach to Compute-and-Forward Network Coding Design\n"+40*"-")
        
    par_dict = {
        "pstart": 0,    # dB/10
        "pend": 2,      # dB/10
        "pnum": 200, # number of points between pstart and pend
        # Lstart, Lend, Lstep only used for standard channels
        "Lstart": 4,    
        "Lend": 4,
        "Lstep": 2,
        # fix channel: np.array([…]); for std gauß channel: anything with size() = 1
        "h": np.array([1.2,0.3,0.8,2.1]),# np.array([None]),  # when None is used -> standard normal distribution will be used to create h's
        "its": 10, # number of iterations per setting (h will be computed newly if std norm is used)
        "time_scale": 1000, # 1000 for ms, 1 for s ...
        "calc_ref": False,
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

















