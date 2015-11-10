# -*- coding: utf-8 -*-


import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import qpr_fundamentials as fun   # fundamential functions like norm2
import qpr_algorithm as alg


def choose_Ku(L):
    """
    Returns a upper bound Ku as listed in Table 1 in [1]
    """
    L_Ku_dict = {2:2, 3:3, 4:4, 5:5, 6:5, 7:5, 8:6, 9:6, 10:6, 11:6, 12:7, 13:6, 14:6, 15:6, 16:6}
    if (L in range(2, 17, 1)):
        return L_Ku_dict[L]
    else:
        return 0

def plot_comprate_vs_P(PdB, compRate, L=None, Ku=None, ITER_NUM=None):
    """
    Plots the compRate vector versus the Power in dB (PdB). Additional informations are the following parameters
    
    Parameters
    ----------
    L: int
        number of senders in system
    Ku: int
        upper bound in qpr calculation
    ITER_NUM: int
        number of iterations/calculation at one power value in qpr
    """
    fig = plt.figure(num=1)
    plt.plot(PdB, compRate, 'go ', label="compRate", figure=fig)
    #plt.plot(Pdb, fmin, 'ro ', label="fmin", figure=fig)
    plt.xlabel("P(dB)")
    plt.ylabel("Average computation rate(bits/channel use)")
    if (L==None) and (Ku==None) and (ITER_NUM==None):
        title_str = ""
    else:
        title_str = "L={}, Ku={}, Iterations per sample={}"
    plt.title(title_str.format(L, Ku, ITER_NUM))
    plt.legend(loc='lower right')

def plot_curve_fit(PdB, ydata_fit):
    """
    Plots the fitted curve into the current figure.
    """
    fig = plt.gcf()
    plt.plot(PdB, ydata_fit, 'r-', label="compRate fit", figure=fig)
    plt.legend(loc='lower right')

    
def write_iteration_to_file(f, ITER_RANGE, compRate_arr, aSquare_arr, h_arr, p_val, p_idx, PdB, Ku, L):
    """
    Writes computation rate, aSquare and channel vector and some other informations into a given file f.
    """
    f.write("CoR:\taSquare\t\t" + "| channel Vector h\tP={:.2f}\tPdB={:.1f}\tKu={}\n".format(p_val, PdB[p_idx], Ku ))
    for idx in ITER_RANGE:
        f.write("{:.4f}:\t".format(compRate_arr[idx] ) )
        # make one flat array with aSquare and h of length 2*L, plus one  coloumn for compRate
        a = np.reshape( [ aSquare_arr[idx], h_arr[idx] ] , 2*L )  
        linestr = L*"{:.0f} " + "\t|" + L*"{:.4f}\t" + "\n"
        f.write( linestr.format( *a ) )
    f.write("\n")

def comp_rate(h, a, P):
    """
    Calculates the Computation Rate as known as the function:
    R(h, a) = 1/2 * log2(1/ (norm2(a) - (P*(h*a)^2))/(1+ P*norm2^2(h)) )
    """
    ha_dotprod_pow2 = np.power(np.dot(h, a), 2)
    denominator = fun.norm2(a) - ( P * ha_dotprod_pow2 / (1 + P * fun.norm2(h) ) )
    R = 0.5 * log2_plus(1 / denominator)
    return R 

def log2_plus(v):
    """
    Returns the log in respect to base 2, if this is bigger or equal 0, whereas the function returns 0 if it is smaller than 0. 
    """
    log2_v = np.log2(v)
    if log2_v < 0:
        return_val = 0
    else:
        return_val = log2_v
    return return_val

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
    plot_comprate_vs_P(PdB, compRate_arr)

def qpr_lin_fit(x, m, n):
    """
    Linear function model to fit the measuremnt data.
    """
    return x*m + n
 
def qpr_main(h=None, L=2, do_plot=False, do_dump=False, do_curve_fit=False, ITER_NUM=10):
    """
    Main function to do quadratic programming relaxation as described in
    [ZM14] 'A quadratic Programming Relaxation Approach Compute-and-Forward Network Coding Design' by
    Zhou and Mow, 2014, Clear Water Bay, Kowloon , Hong Kong
    
    h: np.array   
        channel vector (if None np.random.standard_normal(L) is used)
    L: int  
        length of channel vector (number senders)
    do_plot: bool
        if True, plots are done
    do_dump: bool
        if True, equation coefficients, channel vectors, computation Rate and more infos were dumped in a file
    do_curve_fit: bool
        if True, a curve fitting algorithm is done to calculate a linear curve. If do_plot==True the fitted curve is plotted into the figure 
    ITER_NUM: int
        Number of iterations. Only used if h == None. Then the channel vector is standard_normal distributed
    """
    if do_dump:
        f = open("aSquare_dump.txt", "w")
    if h==None:
        use_standard_normal = True
    else:
        use_standard_normal = False
        ITER_NUM = 1
    ITER_RANGE =range(0,ITER_NUM)# [0]    # only one iteration is necessary, because the result would always be the same
        
    Ku = choose_Ku(L)  # upper bound K
    
    
    compRate_arr = np.zeros(ITER_NUM)  # esult array to hold compRate
    aSquare_arr =  np.zeros([ITER_NUM, L])  # result array to hold aSquare
    h_arr = np.zeros([ITER_NUM, L])  
    
    P_NUM = 50
    P = np.logspace(0, 2, num=P_NUM)
    PdB = 10 * np.log10(P)    # relative to Noise, Noise = 1
    compRate_average = np.zeros(P_NUM)
    
    # iteration loop
    for p_idx, p_val in enumerate(P):
    
        for n in ITER_RANGE:
            if use_standard_normal:
                h_arr[n] = np.absolute( np.random.standard_normal(L) ) # generating a gaussion standard normal distributed
            else:
                h_arr[n] = h
            # run quadtratic programming relaxation algorithm
            # result is a tupel of two arrays
            (compRate_arr[n], aSquare_arr[n]) = alg.qp_relax(h_arr[n], p_val, Ku, L)          
        
        if do_dump:
            write_iteration_to_file(f, ITER_RANGE, compRate_arr, aSquare_arr, h_arr, p_val, p_idx, PdB, Ku, L)
        
        compRate_average[p_idx] = np.average(compRate_arr, axis=0) # average compRate
        
        # evaluation progress print out
        print "{}:\t{:.4f}\t{:.4f}\t{:.4f}\t".format(p_idx, p_val, PdB[p_idx], compRate_average[p_idx])
    
    # post iteration work
    if do_dump:
        f.close()  
    if do_plot:      
        plot_comprate_vs_P(PdB, compRate_average, L, Ku, ITER_NUM)
    if do_curve_fit:
        (fit_parameters, fit_parameters_cov) = curve_fit(qpr_lin_fit, PdB, compRate_average)
        if do_plot:
            m = fit_parameters[0]
            n = fit_parameters[1]
            plot_curve_fit(PdB, qpr_lin_fit(PdB, m, n))

"""
MAIN
"""
if __name__ == '__main__':
    qpr_main(L=16, do_plot=True, do_dump=False, do_curve_fit=True, ITER_NUM=1000)
    #comp_rate_reference(np.array( [0.4,1.3,1.2,0.6] ) , np.array( [1, 4, 4, 2]) ) 
    
    

"""
[Computation rate: Koeffizienten a]

h = [0.4,1.3,1.2,0.6] mit Leistung P = 100
---------------------
1.3044: (1, 2, 2, 1)
1.0447: (1, 4, 4, 2)
0.9729: (1, 3, 3, 1)
0.8871: (0, 1, 1, 0)
0.7248: (2, 6, 6, 3)
0.7088: (0, 1, 1, 1)

h = [1.2,0.3,0.8,2.1] mit Leistung P = 100
---------------------
1.3530: (3, 1, 2, 5)
1.2914: (1, 0, 1, 2)
1.1123: (4, 1, 3, 7)
0.9009: (3, 1, 2, 6)
0.8565: (2, 0, 1, 3)
0.8046: (5, 1, 3, 8)

h = [0.1,0.3,0.8,0.5] mit Leistung P = 100
---------------------
1.0294: (0, 1, 2, 1)
0.8448: (0, 1, 3, 2)
0.8448: (0, 0, 1, 1)
0.7370: (0, 0, 1, 0)
0.5922: (0, 1, 1, 1)
0.4183: (1, 2, 5, 3)
""" 