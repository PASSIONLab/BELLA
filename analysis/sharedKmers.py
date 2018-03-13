import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

#
# Estimation of the nonzero counts of the overlap matrix for non-repetitive genomes
# NOTE: k-mers are found in the same location
#
def func(L, p, k):
	p = pow(p,2)
	S = 0.0
	for i in range(L-k):
		S = S + i*L*pow(1-p,2)*pow(p,k+i)
	return S

print(func(L=700, p=0.85, k=17))

#
# Plot probability function
#
c = [] 	# overlap length array
d = []	# probability array
for i in range(15,8000):
	d.append(func(L=i, p=0.85, k=15))
	c.append(i)

x = [] 	# overlap length array
y = []	# probability array
for i in range(17,8000):
	y.append(func(L=i, p=0.85, k=17))
	x.append(i)

z = [] 	# overlap length array
w = []	# probability array
for i in range(19,8000):
	w.append(func(L=i, p=0.85, k=19))
	z.append(i)

a = [] 	# overlap length array
b = []	# probability array
for i in range(21,8000):
	b.append(func(L=i, p=0.85, k=21))
	a.append(i)


plt.plot(c,d, '#87CEFA', linewidth = 2.0, label='K = 15')
plt.plot(x,y, '#3CB371', linewidth = 2.0, label='K = 17') 
plt.plot(z,w, '#FF7F50', linewidth = 2.0, label='K = 19') 
plt.plot(a,b, '#DC143C', linewidth = 2.0, label='K = 21') 
plt.title("Expected number of shared k-mers S given an overlap length L", size=12.0)
plt.ylabel('S', size=12.0)
plt.xlabel('Overlap length', size=12.0)
plt.show()

#
# End of program
#