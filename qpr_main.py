# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
#import qpr_fundamentials as fun   # fundamential functions like norm2
import qpr_algorithm as alg

def choose_Ku(L): # Table from paper
    L_Ku_dict = {2:2, 3:3, 4:4, 5:5, 6:5, 7:5, 8:6, 9:6, 10:6, 11:6, 12:7, 13:6, 14:6, 15:6, 16:6}
    if (L in range(2, 17, 1)):
        return L_Ku_dict[L]
    else:
        return 0

def plot_comprate_vs_P(PdB, compRate, L, Ku, ITER_NUM):
    fig = plt.figure(num=1)
    plt.plot(PdB, compRate, 'go ', label="compRate", figure=fig)
    #plt.plot(Pdb, fmin, 'ro ', label="fmin", figure=fig)
    plt.xlabel("P(dB)")
    plt.ylabel("Average computation rate(bits/channel use)")
    plt.title("L={}, Ku={}, Iterations per sample={}".format(L, Ku, ITER_NUM))
    plt.legend()

def qpr_main():
    
    L = 8  # number of channels   
    Ku = choose_Ku(L)  # upper bound K
    
    ITER_NUM = 10
    ITER_RANGE = range(0,ITER_NUM)
    compRate_arr = np.zeros(ITER_NUM)  # esult array to hold compRate
    aSquare_arr =  np.zeros([ITER_NUM, L])  # result array to hold aSquare
    h_arr = np.zeros([ITER_NUM, L])  
    
    P_NUM = 50
    P = np.logspace(0, 2, num=P_NUM)
    PdB = 10 * np.log10(P)    # relative to Noise, Noise = 1
    compRate_average = np.zeros(P_NUM)
    
    f = open("aSquare_dump.txt", "w")
    
    for p_idx, p_val in enumerate(P):
    
        for n in ITER_RANGE:
    
            h_arr[n] = np.absolute( np.random.standard_normal(L) ) # generating a gaussion standard normal distributed 

            # run quadtratic programming relaxation algorithm
            # result is a tupel of two arrays
            res_tupel = alg.qp_relax(h_arr[n], p_val, Ku, L) 
            # allocate result tupel            
            compRate_arr[n] = res_tupel[0]  # array containing [compRate, fmin]
            aSquare_arr[n] = res_tupel[1]   # aSquare
        
        f.write("aSquare\tchannelVector\n")
        for idx in ITER_RANGE:
            f.write("%s\n" % [aSquare_arr[idx], h_arr[idx]] )
        
        
        compRate_average[p_idx] = np.average(compRate_arr, axis=0) # calculate average compRate and fmin

        print "{}:\t{:.4f}\t{:.4f}\t{:.4f}\t".format(p_idx, p_val, PdB[p_idx], compRate_average[p_idx])
    plot_comprate_vs_P(PdB, compRate_average, L, Ku, ITER_NUM)
    

'''
MAIN
'''
if __name__ == '__main__':
    qpr_main()
    