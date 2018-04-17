#
# @author: Giulia Guidi
#

import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

#
# Markov Chain transition matrix to get the probability of having (k) consecutive correct bases on two different 
# sequences given an overlap (t) between the two sequences and the probability of having a correct base (p)
# NOTE: k-mers are found in the same location
#
def func(overlap, p, k):
    transixMat = np.zeros((k+1, k+1))
    for i in range(k):
        transixMat[i, 0] = (1-p*p)
        transixMat[i, i + 1] = p*p
    transixMat[k, k] = 1.0
    return np.linalg.matrix_power(transixMat, overlap)[0, k]

print(func(overlap=500, p=0.85, k=17))

#
# Plot probability function
#
c = [] 	# overlap length array
d = []	# probability array
for i in range(15,10000):
	d.append(func(overlap=i, p=0.85, k=15))
	c.append(i)

x = [] 	# overlap length array
y = []	# probability array
for i in range(17,10000):
	y.append(func(overlap=i, p=0.85, k=17))
	x.append(i)

z = [] 	# overlap length array
w = []	# probability array
for i in range(19,10000):
	w.append(func(overlap=i, p=0.85, k=19))
	z.append(i)

a = [] 	# overlap length array
b = []	# probability array
for i in range(21,10000):
	b.append(func(overlap=i, p=0.85, k=21))
	a.append(i)


plt.plot(c,d, '#87CEFA', linewidth = 2.0, label='K = 15')
plt.plot(x,y, '#3CB371', linewidth = 2.0, label='K = 17') 
plt.plot(z,w, '#FF7F50', linewidth = 2.0, label='K = 19') 
plt.plot(a,b, '#DC143C', linewidth = 2.0, label='K = 21') 
plt.title("Probability to find at least one shared k-mer given the overlap length L", size=12.0)
plt.ylabel('P(L)', size=12.0)
plt.xlabel('Overlap length', size=12.0)
plt.show()

#
# End of program
#