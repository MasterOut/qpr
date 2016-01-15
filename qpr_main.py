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
    filename = 'qpr_run_' + date_str + '.txt'
#    fieldnames = ['nr', 'P', 'Pdb', 'h', 'a', 'a_occ', 'R', 'time', 'R_ref']
    w = qcsv.csv_dict_writer(filename)    
    
    # generate P and PdB
    P = np.logspace(par_dict["pstart"] /10, par_dict["pend"] /10, par_dict["pnum"])
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
                
            if par_dict["h_absolute"]:
                    h = np.absolute(h) 
                
            cr = np.zeros((its,1))  # array to hold cr values
            a = np.zeros((its, L))  # array to hold a coefficients
    
            """ iteration loop with time measurement"""
            tstart = time.time()
            for i in its_range:
                cr[i], a[i] =  alg.qp_relax(h[i], p, Ku, L)
                #if np.dot(h[i], a[i]) < 0:
                #    print "h[i] dot a[i] < 0 !"
            tend = time.time()
            time_needed = (tend - tstart) * time_scale
            
            # cnt appearance for each choefficient vector
            a_list, a_occ_list = fun.cnt_appearance(a)
            cr_mean = np.mean(cr)
            cr_std = np.std(cr)
            #cr_max = np.max(cr)
            #cr_min = np.min(cr)
            
            """ Comparison with reference computation rate formula """
            if par_dict["calc_ref"]:
                a_most = a_list[0] # most calculated coeff vector
                for i_h, v_h in enumerate(h):
                    cr_ref[i_h] = fun.comp_rate(v_h, a_most, p)
                cr_ref_av = np.average(cr_ref)
            #print "\t{0}\t{1}\t{2}\t{3}".format(cr_mean, cr_max, cr_min, cr_ref_av)
            
            """ Prepare dump dictionary """
            w_dict = {'nr': ind,'P': p,'Pdb': Pdb[ind], 'a': a_list,'a_occ': a_occ_list, 'Rmean': cr_mean, 'Rstd': cr_std, 'time': time_needed,'R_ref_av': cr_ref_av}
            if use_std_normal:
                w_dict['h'] = 'std'
            else:
                w_dict['h'] = h
            """ Write dictionary to file """
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
        "pstart": 0,    # dB
        "pend": 20,     # dB
        "pnum": 100,    # number of points between pstart and pend
        # Lstart, Lend, Lstep only used for standard channels
        "Lstart": 2,    # channel and coefficient length START
        "Lend": 16,     # END
        "Lstep": 2,     # Stepwidth
        "h": np.array([None]),  # None for gaussian channel vector
        "h_absolute": True,     # True: only positive values in channel vector
        "its": 1000,     # number of iterations per setting
        "time_scale": 1000, # 1000 for ms, 1 for s ...
        "calc_ref": False,  # 
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

















