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

def func(L, p, k):
    states = k+1
    transixMat = np.zeros((states, states))
    for i in range(k):
        transixMat[i, 0] = (1-p*p)
        transixMat[i, i + 1] = p*p
    transixMat[k, k] = 0.0

    v = np.zeros(states)
    v[0] = 1
    vk   = []
    for s in range(L):
        v = np.dot(v, transixMat)
        vk.append(v[k])
    # print(vk)

    transixMat[k, k] = 1.0
    wk = 0.0
    for s in range(L):
        w = np.zeros(states)
        w[0] = vk[L-s-1]
        w = np.dot(w, np.linalg.matrix_power(transixMat, s))
        wk = wk + w[k]
    return np.linalg.matrix_power(transixMat, L)[0, k], wk

ksize=17
pcorr=0.85

# p1, p2 = func(L=2000, p=pcorr, k=ksize)
# print(p1, p2)

probs = [func(L=i, p=pcorr, k=ksize) for i in range(500, 15005, 500)]
probs = np.asarray(probs)
# print(probs)
# print(probs[:,0], probs[:,1])

plt.plot(np.arange(500, 15005, 500), probs[:,0], color='black', label='BMC')
plt.plot(np.arange(500, 15005, 500), probs[:,1], color='red',   label='2-states NMC')
# plt.plot(np.arange(1000, 30005, 1000), probs[:,0]**2, color='blue', label='BMC**2')
plt.title('Nested Markov Chain')
plt.xlabel('Overlap Length')
plt.ylabel('Cumulative Probability')
plt.legend()
plt.show()

#
# End of program
#