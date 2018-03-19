import numpy as np
# import matplotlib.patches as mpatches
# import matplotlib.pyplot as plt

#
# Markov Chain transition matrix to get the probability of having (k) consecutive correct bases on two different 
# sequences given an overlap (t) between the two sequences and the probability of having a correct base (p)
# NOTE: k-mers can be found in different locations
#
# The following two functions allow to move from multi-dimensional matrix to bidimensional matrix:
# 1) Function to pass from the matrix index to the actual state on sequence (i)
def idxtoSuccessOn1(idx,k):
	return int(idx/(k+1))
# 2) Function to pass from the matrix index to the actual state on sequence (j)
def idxtoSuccessOn2(idx,k):
	return idx%(k+1)

def statetoIdx(i,j,k):
	return i*(k+1)+j
#
# Info needed during the analyitical computation of the general form
#
# class matInfo:
# 	# identificative tuple
#     name = () 
#     # entry value
#     value = 0.0 

#
# Function to build the transition matrix according to Markov Chain Model 
#
def func(t,p,k):
	nstates = pow(k+1,2)-1; # number of states of the transition matrix
	transixMat = np.zeros((nstates+1, nstates+1)) # initialization to zero
	#
	# The rows of this two-dimensional matrix represent the starting states and cols the ending states, each cell represents the probability of moving from one state to another
	# Each index represent a state which encloses the condition of both the sequences (number of successes for one sequence (i) and for the other one (j))
	#
	for start in range(nstates+1):
		for end in range(nstates+1):
			# Do not consider the cases where the row or the col is equal to K as those cases have different formulas
			if idxtoSuccessOn1(start,k) != k and idxtoSuccessOn2(start,k) != k:
				# transixMat[0,0]
				if end == 0: 
					# Probability that both the bases are incorrectly sequenced
					transixMat[start,end] = (1-p)*(1-p) 
				# transixMat[i+1][0]
				if idxtoSuccessOn1(end,k) == idxtoSuccessOn1(start,k)+1 and idxtoSuccessOn2(end,k) == 0:
					# Probability to correctly sequence a base just on the (i) sequence
					transixMat[start,end] = p*(1-p)
				# transixMat[0][j+1]
				if idxtoSuccessOn1(end,k) == 0 and idxtoSuccessOn2(end,k) == idxtoSuccessOn2(start,k)+1:
					# Probability to correctly sequence a base just on the (j) sequence
					transixMat[start,end] = (1-p)*p
				# transixMat[i+1][j+1] with i != k and j != k
				if idxtoSuccessOn1(end,k) == idxtoSuccessOn1(start,k)+1 and idxtoSuccessOn2(end,k) == idxtoSuccessOn2(start,k)+1:
					# Probability to correctly sequence a base on both the sequences
					transixMat[start,end] = p*p
			# Cases where sequence (j) got k successes
			elif idxtoSuccessOn1(start,k) != k and idxtoSuccessOn2(start,k) == k:
				# transixMat[0][k]
				if idxtoSuccessOn1(end,k) == 0 and idxtoSuccessOn2(end,k) == k:
					# (j) sequence is already arrived at k successes, so successes cannot go back to zero
					# The base on (i) sequence can be incorrectly sequenced with probability (1-p) and send back the successes to zero
					transixMat[start,end] = 1-p
				# transixMat[i+1][k]
				if idxtoSuccessOn1(end,k) == idxtoSuccessOn1(start,k)+1 and idxtoSuccessOn2(end,k) == k:
					# (j) sequence is already arrived at k successes, so it cannot achieve other successes
					# The base on (i) sequence can be correctly sequenced with probability p and increase the successes
					transixMat[start,end] = p 
			# Cases where sequence (i) got k successes
			elif idxtoSuccessOn1(start,k) == k and idxtoSuccessOn2(start,k) != k:
				# transixMat[k][0]
				if idxtoSuccessOn2(end,k) == 0 and idxtoSuccessOn1(end,k) == k:
					# (i) sequence is already arrived at k successes, so successes cannot go back to zero
					# The base on (j) sequence can be incorrectly sequenced with probability (1-p) and send back the successes to zero
					transixMat[start,end] = 1-p
				# transixMat[k][j+1]
				if idxtoSuccessOn2(end,k) == idxtoSuccessOn2(start,k)+1 and idxtoSuccessOn1(end,k) == k:
					# (i) sequence is already arrived at k successes, so it cannot achieve other successes
					# The base on (j) sequence can be correctly sequenced with probability p and increase the successes
					transixMat[start,end] = p				

	# Both sequences obtained k successes
	transixMat[nstates,nstates] = 1
	# Exponantiation 
	finalMat = np.linalg.matrix_power(transixMat, t)
	# print(t)

	#	if 0 in finalMat:
	#		print("There is at least a 0, so it's not a regular Markov chain")
	#	
	#	globalList = [] # global
	#	print(general(finalMat,nstates,t,globalList))

	return np.linalg.matrix_power(transixMat, t)[0, nstates]

#
# A^t[0,(k+1)^2] = SUM from i = 0 to (k+1)^2 A^(t-1)[0,i]*A^(t-1)[i,(k+1)^2-1]
# if t = 1, then set all products with the corresponding A = 0 to 0
#
# def general(transixMat,nstates,t,globalList):
# 	localInfo = matInfo()
# 	localInfo.value = 0.0 # back to zero at each function call
# 	print(t)
# 	if t == 1:
# 		# set all products with the corresponding A = 0 to 0
# 		for item in globalList:
# 			i = item.name[1]
# 			j = item.name[2]
# 			if transixMat[i,j] == 0.0:
# 				item.value = 0.0
# 		return tr
# 	else:
# 		for i in range(nstates+1):
# 			for j in range(nstates+1):
# 				localInfo.name = (t,i,j)
# 				for k in range(nstates+1):
# 					# print(str(type(general(transixMat,nstates,t-1,globalList))))
# 					localInfo.value += general(transixMat,nstates,t-1,globalList)[i, k]*general(transixMat,nstates,t-1,globalList)[k, j]
# 					# np.linalg.matrix_power(transixMat, t)[0, nstates] += localInfo.value
# 				print(localInfo.value)
# 		globalList.append(localInfo)
# 		return np.linalg.matrix_power(transixMat, t)

# t = overlap length
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