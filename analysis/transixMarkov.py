import numpy as np
# import matplotlib.patches as mpatches
# import matplotlib.pyplot as plt

#
# Markov Chain transition matrix to get the probability of having (k) consecutive correct bases on two different 
# sequences given an overlap (t) between the two sequences and the probability of having a correct base (p)
#
# The following two functions allow to move from multi-dimensional matrix to bidimensional matrix:
# 1) Function to pass from the matrix index to the actual state on sequence (i)
def idxtoRow(idx,k):
	return int(idx/(k+1))
# 2) Function to pass from the matrix index to the actual state on sequence (j)
def idxtoCol(idx,k):
	return idx%(k+1)
#
# Function to build the transition matrix according to Markov Chain Model 
#
def func(t, p, k):
	nstates = pow(k+1,2)-1; # number of states of the transition matrix
	transixMat = np.zeros((nstates+1, nstates+1)) # initialization to zero
	#
	# The rows of this two-dimensional matrix represent the starting states and cols the ending states, each cell represents the probability of moving from one state to another
	# Each index represent a state which encloses the condition of both the sequences (number of successes for one sequence (i) and for the other one (j))
	#
	for start in range(nstates+1):
		for end in range(nstates+1):
			# Do not consider the cases where the row or the col is equal to K as those cases have different formulas
			if idxtoRow(start,k) != k and idxtoCol(start,k) != k:
				# transixMat[0,0]
				if end == 0: 
					# Probability that both the bases are incorrectly sequenced
					transixMat[start,end] = (1-p)*(1-p) 
				# transixMat[i+1][0]
				if idxtoRow(end,k) == idxtoRow(start,k)+1 and idxtoCol(end,k) == 0:
					# Probability to correctly sequence a base just on the (i) sequence
					transixMat[start,end] = p*(1-p)
				# transixMat[0][j+1]
				if idxtoRow(end,k) == 0 and idxtoCol(end,k) == idxtoCol(start,k)+1:
					# Probability to correctly sequence a base just on the (j) sequence
					transixMat[start,end] = (1-p)*p
				# transixMat[i+1][j+1] with i != k and j != k
				if idxtoRow(end,k) == idxtoRow(start,k)+1 and idxtoCol(end,k) == idxtoCol(start,k)+1:
					# Probability to correctly sequence a base on both the sequences
					transixMat[start,end] = p*p
			# Cases where sequence (j) got k successes
			elif idxtoRow(start,k) != k and idxtoCol(start,k) == k:
				# transixMat[0][k]
				if idxtoRow(end,k) == 0 and idxtoCol(end,k) == k:
					# (j) sequence is already arrived at k successes, so successes cannot go back to zero
					# The base on (i) sequence can be incorrectly sequenced with probability (1-p) and send back the successes to zero
					transixMat[start,end] = 1-p
				# transixMat[i+1][k]
				if idxtoRow(end,k) == idxtoRow(start,k)+1 and idxtoCol(end,k) == k:
					# (j) sequence is already arrived at k successes, so it cannot achieve other successes
					# The base on (i) sequence can be correctly sequenced with probability p and increase the successes
					transixMat[start,end] = p 
			# Cases where sequence (i) got k successes
			elif idxtoRow(start,k) == k and idxtoCol(start,k) != k:
				# transixMat[k][0]
				if idxtoCol(end,k) == 0 and idxtoRow(end,k) == k:
					# (i) sequence is already arrived at k successes, so successes cannot go back to zero
					# The base on (j) sequence can be incorrectly sequenced with probability (1-p) and send back the successes to zero
					transixMat[start,end] = 1-p
				# transixMat[k][j+1]
				if idxtoCol(end,k) == idxtoCol(start,k)+1 and idxtoRow(end,k) == k:
					# (i) sequence is already arrived at k successes, so it cannot achieve other successes
					# The base on (j) sequence can be correctly sequenced with probability p and increase the successes
					transixMat[start,end] = p				

	# Both sequences obtained k successes
	transixMat[nstates,nstates] = 1
	# Exponantiation 
	return np.linalg.matrix_power(transixMat, t)[0, nstates]

print(func(t=700, p=0.85, k=17))

#
# Plot probability function
#
# c = [] 	# overlap length array
# d = []	# probability array
# for i in range(15,700):
# 	d.append(func(t=i, p=0.85, k=15))
# 	c.append(i)
#
# x = [] 	# overlap length array
# y = []	# probability array
# for i in range(17,700):
# 	y.append(func(t=i, p=0.85, k=17))
# 	x.append(i)
# 
# z = [] 	# overlap length array
# w = []	# probability array
# for i in range(19,700):
# 	w.append(func(t=i, p=0.85, k=19))
# 	z.append(i)
# 
# a = [] 	# overlap length array
# b = []	# probability array
# for i in range(21,700):
# 	b.append(func(t=i, p=0.85, k=21))
# 	a.append(i)
# 
# 
# plt.plot(c,d, '#87CEFA', linewidth = 2.0, label='K = 15')
# plt.plot(x,y, '#3CB371', linewidth = 2.0, label='K = 17') 
# plt.plot(z,w, '#FF7F50', linewidth = 2.0, label='K = 19') 
# plt.plot(a,b, '#DC143C', linewidth = 2.0, label='K = 21') 
# plt.title("Probability of having one correct k-mer on both the sequences (K) given an overlap region", size=12.0)
# plt.ylabel('P(K)', size=12.0)
# plt.xlabel('Overlap length', size=12.0)
# plt.show()

#
# End of program
#