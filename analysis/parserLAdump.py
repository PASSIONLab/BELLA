#!/usr/bin/python
#
# Reformats the output of 'LAdump -cl' (provided in a file specified via argument 1)
# in bella format with one exception: read ID's are output in place of read names.
# The unique set of read IDs are recorded and output for lookup in the corresponding DAZZ_DB database instance.
# See translateDalignerOut.sh for the full translation procedure (which includes replacing the ID's with names).
#
# @author:mme
#
from sys import argv

if len(argv) < 2:
  print("1 argument is required: a file containing 'LAdump -cl' formatted alignment output")
  quit()

ofilename=argv[1]+".bella"
idset=set()
with open(argv[1],'r') as infile:
  line=[]
  for l in infile:
    lsplit = l.split()
    if (lsplit[0]=='P'):
      line+=lsplit[1:4]
      idset=idset | set(lsplit[1:3]) #union
    if (lsplit[0]=='L'):
      line+=lsplit[1:]
    if (lsplit[0]=='C'):
      # r1 r2 c/n [3] [4] r1L [6] [7] r2L
      line.insert(3, lsplit[1])
      line.insert(4, lsplit[2])
      line.insert(6, lsplit[3])
      line.insert(7, lsplit[4])
      with open(ofilename, 'a') as ofile:
        lineout=' '.join(x for x in line)
        ofile.write(lineout+'\n')
      line=[] # reset to process next 3 lines

# write the set of ids, 1 per line, to a file to use with DBshow
idsfilename=argv[1]+".indices.txt"
with open(idsfilename,'w') as ofile:
  for id in idset:
    ofile.write(id+'\n')
