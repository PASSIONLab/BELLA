#
# @author: Giulia Guidi
#

import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

#
# Nested Markov Chain transition matrix to get the probability of having 2 (k) consecutive correct k-mers on two different 
# sequences given an overlap (t) between the two sequences and the probability of having a correct base (p)
# NOTE: k-mers are found in the same location
#

def nmc(L, p, k, w):
    states = w*k+1
    transixMat = np.zeros((states, states))
    for i in range(w):
        for j in range(k):
            transixMat[j+i*k, i*k]       = (1-p*p)
            transixMat[j+i*k, j+i*k + 1] = p*p      
    transixMat[k*w, k*w] = 1.0
    # print(transixMat)
    return np.linalg.matrix_power(transixMat, L)[0, w*k]

ksize=17
pcorr=0.85

## data generation

# probs1 = [nmc(L=i, p=pcorr, k=ksize, w=1) for i in range(2000, 13005, 500)]
# probs1 = np.asarray(probs1)

# probs2 = [nmc(L=i, p=pcorr, k=ksize, w=2) for i in range(2000, 13005, 500)]
# probs2 = np.asarray(probs2)

# probs3 = [nmc(L=i, p=pcorr, k=ksize, w=3) for i in range(2000, 13005, 500)]
# probs3 = np.asarray(probs3)

# probs4 = [nmc(L=i, p=pcorr, k=ksize, w=4) for i in range(2000, 13005, 500)]
# probs4 = np.asarray(probs4)

# probs5 = [nmc(L=i, p=pcorr, k=ksize, w=5) for i in range(2000, 13005, 500)]
# probs5 = np.asarray(probs5)

# probs6 = [nmc(L=i, p=pcorr, k=ksize, w=6) for i in range(2000, 13005, 500)]
# probs6 = np.asarray(probs6)

# plt.hlines(0.90, 2000, 13000, colors='red', linestyles='--', label='0.90')

# plt.plot(np.arange(2000, 13005, 500), probs1, color='black',        label='BMC')
# plt.plot(np.arange(2000, 13005, 500), probs2, color='royalblue',    label='2-states NMC')
# plt.plot(np.arange(2000, 13005, 500), probs3, color='darkorchid',   label='3-states NMC')
# plt.plot(np.arange(2000, 13005, 500), probs4, color='green',        label='4-states NMC')
# plt.plot(np.arange(2000, 13005, 500), probs5, color='darkorange',   label='5-states NMC')
# plt.plot(np.arange(2000, 13005, 500), probs6, color='firebrick',    label='6-states NMC')
# plt.title('Nested Markov Chain')
# plt.xlabel('Overlap Length')
# plt.ylabel('Cumulative Probability')
# plt.legend()
# plt.show()

#
# End of program
#