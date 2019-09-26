#
# @author: Giulia Guidi
#

import numpy as np

#
# Markov Chain transition matrix to get the probability of having (k) consecutive correct bases on two different 
# sequences given an overlap (t) between the two sequences and the probability of having a correct base (p)
# NOTE: k-mers are found in the same location
#

#
# Expected value of out Markov chain
#

# Adding the entries in the top row, we get the expected number of steps
def getResult(fundamentalMat):
	return fundamentalMat[1, :].sum()

# The matrix N is the fundamental matrix for P
def genN(transixMatQ, k):
	I = np.identity(k)
	tmp = I - transixMatQ
	fundamentalMat = np.linalg.inv(tmp)
	return fundamentalMat

# Let Q be the sub-matrix of P without the rows and columns of any absorbing states
def genQ(transixMat, k):
	transixMatQ = transixMat[:-1, :-1]
	return transixMatQ

def myOverlap(p, k):
	transixMat = np.zeros((k+1, k+1))
	for i in range(k):
		transixMat[i, 0] = (1-p*p)
		transixMat[i, i + 1] = p*p
	transixMat[k, k] = 1.0
	transixMatQ = genQ(transixMat, k)
	fundamentalMat = genN(transixMatQ, k)
	return getResult(fundamentalMat)

#
# End of program
#