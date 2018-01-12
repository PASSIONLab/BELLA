import csv 
from sets import Set
import sys, os, fileinput, string

with open('BELLAreads.txt', 'r') as bellafile:
	bellareader = csv.reader(bellafile, delimiter=',')
	BELLApair = []
	pair = []
	for line in bellareader:
		pair.append(line[0])
		pair.append(line[1])
		BELLApair.append(pair)
		del pair[:]

BWApair = []
out_file = open('BWAMEMoverlap.txt', 'w')

with open(sys.argv[1], 'r') as inputfile: # take BWAMEMread.txt as input 
	treader = csv.reader(inputfile, delimiter=',')
	for tline in treader:
		tstart = int(tline[1])
		tend = int(tline[2])
		inputfile.seek(0,0)
		qreader = csv.reader(inputfile, delimiter=',')
		for qline in qreader:
			qstart = int(qline[1])
			qend = int(qline[2])
			if tstart < qstart:
				if tend < qstart:
					aln = min(tend-qstart, qend-qstart)
			elif tstart > qstart:
				if qend > tstart:
					aln = min(qend-tstart, tend-tstart)
			else:
				aln = min(qend-tstart, tend-tstart)
			if aln > 2000:
				pair = []
				pair.append('@' + tline[0])
				pair.append('@' + qline[0])
				BWApair.append(pair)
				del pair[:]

with open(sys.argv[2], 'r') as inputfilerv: # take BWAMEMread.txt as input 
	treader = csv.reader(inputfilerv, delimiter=',')
	for tline in treader:
		tstart = int(tline[1])
		tend = int(tline[2])
		inputfilerv.seek(0,0)
		qreader = csv.reader(inputfilerv, delimiter=',')
		for qline in qreader:
			qstart = int(qline[1])
			qend = int(qline[2])
			if tstart < qstart:
				if tend < qstart:
					aln = min(tend-qstart, qend-qstart)
			elif tstart > qstart:
				if qend > tstart:
					aln = min(qend-tstart, tend-tstart)
			else:
				aln = min(qend-tstart, tend-tstart)
			if aln > 2000:
				pair = []
				pair.append('@' + tline[0])
				pair.append('@' + qline[0])
				BWApair.append(pair)
				del pair[:]

with open(sys.argv[1], 'r') as inputfile:
	with open(sys.argv[2], 'r') as inputfilerv:
		treader = csv.reader(inputfile, delimiter=',')
		for tline in treader:
			tstart = int(tline[1])
			tend = int(tline[2])
			inputfilerv.seek(0,0)
			qreader = csv.reader(inputfilerv, delimiter=',')
			for qline in qreader:
				qstart = int(qline[1])
				qend = int(qline[2])
				if tstart < qstart:
					if tend < qstart:
						aln = min(tend-qstart, qend-qstart)
				elif tstart > qstart:
					if qend > tstart:
						aln = min(qend-tstart, tend-tstart)
				else:
					aln = min(qend-tstart, tend-tstart)
				if aln > 2000:
					pair = []
					pair.append('@' + tline[0])
					pair.append('@' + qline[0])
					BWApair.append(pair)
					del pair[:]

BELLAtuple = [tuple(lst) for lst in BELLApair]
BWAtuple = [tuple(lst) for lst in BWApair]

out_file.close()
inputfile.close()
inputfilerv.close()

intersection = set(BELLAtuple) & set(BWAtuple)
print len(BELLAtuple)	    # reads found by BELLA
print len(BWAtuple) 		# reads found by BWAMEM
print len(intersection)


