#!/usr/bin/python
#
# Convenience script for viewing alignments above a given THRESHOLD (optional argument 2, default THRESHOLD=2000).
# Expects the alignments to be provided (argument 1) in bella format.
# Does NOT change the original input file.  
# Outputs the filtered alignments into a separate file, with the naming convention <input file name>_filtered<THRESHOLD>.
#
# @author:mme
#
from sys import argv

if (len(argv) < 2):
  print('Exactly 1 input file (in bella alignment format) required for filtering. Optionally provide integer threshold as arg 2. Exiting.')
  exit()

THRESHOLD=int(2000)
if (len(argv) > 2):
  THRESHOLD=int(argv[2])

# Constants for indexing a bella/dibella alignment output line
TRUNC=4 #the line index at (inclusive) which to begin reading values (read names, similarity, and score, are not used in overlap length filtering)
R1_ALGN_START=4-TRUNC
R1_ALGN_END=5-TRUNC
R1_READ_END=6-TRUNC
R2_ALGN_START=7-TRUNC
R2_ALGN_END=8-TRUNC
R2_READ_END=9-TRUNC
# output format after:
#str << r1 << " " << r2
#    << data.similarity << " " << data.score << " "
#    << data.r1_al_start << " " <<  data.r1_al_end << " " << data.length_r1 << " "
#    << data.r2_al_start << " " <<  data.r2_al_end << " " << data.length_r2
#
ofilename=argv[1]+"_filtered"+str(THRESHOLD)
with open(argv[1],'r') as alignments:
  for line in alignments:
    vals=[int(x) for x in line.split()[TRUNC:]]
    r1_algn_length=vals[R1_ALGN_END]-vals[R1_ALGN_START]
    r2_algn_length=vals[R2_ALGN_END]-vals[R2_ALGN_START]
    head=min(vals[R1_ALGN_START], vals[R2_ALGN_START])
    tail=min( (vals[R1_READ_END] - vals[R1_ALGN_END]), (vals[R2_READ_END] - vals[R2_ALGN_END]) )
    overlap = head + tail + int((r1_algn_length+r2_algn_length)/2)
    if (overlap >= THRESHOLD):
      with open(ofilename,'a') as ofile:
        ofile.write(line)
