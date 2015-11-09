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

def plot():
    fig = plt.figure(num=1)
    plt.plot(Pdb, compRate, 'go ', Pdb, fmin, 'ro ', figure=fig)
    plt.xlabel("P(dB)")
    plt.ylabel("Average computation rate(bits/channel use)")
    plt.title("L={}, Ku={}, Iterations per sample={}".format(L, Ku, ITER_NUM))

L = 4  # number of channels   
Ku = choose_Ku(L)  # upper bound K

ITER_NUM = 1000
res = np.zeros([ITER_NUM, 2])

P_NUM = 100
P = np.logspace(-1, 3, num=P_NUM)
Pdb = 10*np.log10(P*1000)

compRate = np.zeros(P_NUM)
fmin = np.zeros(P_NUM)

d = 0.0 # global variabel, needed in qp_quantization and fmin

# FILE
f = open("chcoeff_dump.txt", "w")

for p_idx, p_val in enumerate(P):

    for n in range(0,ITER_NUM,1):

        h = np.random.standard_normal(L)
        
        res[n] = alg.qp_relax(h, p_val, Ku, L, d)
    
    res_average = np.average(res, axis=0) # 
    compRate[p_idx] = res_average[0]
    fmin[p_idx] = res_average[1]
    print "{}:\t{:.4f}\t{:.4f}\t{:.4f}\t".format(p_idx, p_val, 10*np.log10(p_val*1000), compRate[p_idx])

plot()
