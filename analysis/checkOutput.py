#=======================================================================
# Title:  Python script to compare diBELLA and BELLA (alignment) output
# Author: G. Guidi
# Date:   11 Feb 2019
#=======================================================================

# Run: python3 checkOutput.py <file1> <file2>
# file1 is supposed to be the query and file2 the reference

import csv
import sys
import os
import random

cwd = os.getcwd()
pair = {}
list1 = []	# Query pair names
list2 = []	# Reference pair names

with open(sys.argv[1], 'r') as one, open(sys.argv[2], 'r') as two:
	# Here, we compare only read names, not alignment scores, etc.
	file1 = csv.reader(one, delimiter='\t')
	for row in file1:
		pair = [(row[0], row[1])]
		list1.append(pair)
		pair = [(row[1], row[0])]	# Check reversed pair
		list1.append(pair)
	file2 	= csv.reader(two, delimiter='\t')
	for row in file2:
		pair = [(row[2], row[3])]	# diBELLA has names in position 2,3
		list2.append(pair)
		pair = [(row[3], row[2])]	# Check reversed pair
		list2.append(pair)

m = 0
e = 0
lines = 0
randid = random.randint(1,100000)

with open(cwd + '/missingPairj' + str(randid) + '.out', 'w') as missingPair:
	for pair in list2:			# BELLA/reference
		lines = lines + 1
		if pair not in list1:	# diBELLA/query
			for a, b in pair:
				missingPair.write(a + '\t' + b + '\n')
				m = m + 1	# number of missing lines

with open(cwd + '/extraPairj' + str(randid) + '.out', 'w') as extraPair:
	for pair in list1:			# diBELLA/query
		if pair not in list2:	# BELLA/reference
			for a, b in pair:
				extraPair.write(a + '\t' + b + '\n')
				e = e + 1	# number of missing lines

print(lines, "lines in BELLA (reference)")
print(m, "missing lines in diBELLA (query)")
print(e, "extra lines in diBELLA (query)")

