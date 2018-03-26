import sys, os, fileinput, string
import csv

print("\nBella's output post-processing to throw out entries with estimated overlap smaller than a threshold")
print("Usage: python3 ovlParser.py <bella-standard-output> <name-output-file> <threshold>")
#print("Start parsing Bella's output")

def func(s1,e1,l1,s2,e2,l2,t):
	diff1 = int(e1)-int(s1)
	diff2 = int(e2)-int(s2)
	left = min(int(s1),int(s2))
	right = min(int(l1)-int(e1),int(l2)-int(e2))
	result = left+right+(diff1+diff1)/2
	if result < t:
		#print("Pair discharged")
		return False
	else:
		#print("Pair kept")
		return True

t = int(sys.argv[3])
print("Threshold: ", t)
print("Input: ", str(sys.argv[1]))
print("Output: ", str(sys.argv[2]), "\n")

fout = open(sys.argv[2], 'w')
writer = csv.writer(fout, delimiter='\t')
with open(sys.argv[1], 'r') as f:
	reader = csv.reader(f,delimiter='\t')
	for row in reader:
		I1 = row[0] # read name 1
		I2 = row[1]	# read name 2
		SK = row[2]	# num shared kmer
		AS = row[3]	# alignment score
		S1 = row[4]	# start read 1
		E1 = row[5]	# end read 1
		L1 = row[6]	# length read 1
		S2 = row[7]	# start read 2
		E2 = row[8]	# end read 2
		L2 = row[9]	# length read 2
		if(func(s1=S1,e1=E1,l1=L1,s2=S2,e2=E2,l2=L2,t=t)):
			writer.writerow(row)
fout.close()
#print("End parsing Bella's output")