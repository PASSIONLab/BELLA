import sys
import os
import numpy as np
import scipy.special
import matplotlib.pyplot as plt

k = int(sys.argv[1])
d = int(sys.argv[3])
e = float(sys.argv[2])

def myprob(m,d,k,e):
	bincoef = scipy.special.comb(d, m, exact=True)

	a = np.power(1-e, k*m)
	b = np.power(1-e, k)
	c = np.power(1-b, d-m)

	probability = bincoef * a * b * c
	return probability

x = []
y = []
for m in range(d):
	x.append(m)
	y.append(myprob(m,d,k,e))


plt.bar(x,y)
plt.ylabel('Probability of Correct K-mer')
plt.xlabel('K-mer Frequency')
plt.title('Correct K-mer Probability - ' + str(k) + '-' + str(e) + '-' + str(d))

plt.tight_layout()
plt.savefig('kmer-probability-' + str(k) + '-' + str(e) + '-' + str(d) + '.pdf')

plt.show()

