#!/usr/bin/python
from sys import argv
from fileinput import FileInput

if len(argv) < 4:
  print("3 arguments are required: (1) a text file containing a list of reads by name, 1 per line, (2) a text file containing a list of ID's for (1) in the corresponding order, and (3) a file containing the corresponding BELLA/diBELLA untagged alignment output")
  quit()

tag_map=dict()
namesfile=argv[1]
idsfile=argv[2]
with open(idsfile) as ids:
  with open(namesfile) as names:
     tag_map={k.rstrip():(v.split()[:1])[0].rstrip()[1:] for (k,v) in zip(ids,names)}

#print(tag_map)

#open alignments file and replace first two words with of each line with value from dictionary
space=' '
filename=argv[3]
with FileInput(filename, inplace=True, backup='.bak') as file:
  for line in file:
    linesplit=line.split()
    linesplit[0]=tag_map[linesplit[0]].rstrip()
    linesplit[1]=tag_map[linesplit[1]].rstrip()
    lineout=' '.join(x for x in linesplit)
    print(lineout)
