from simplesam import Reader, Writer
import inspect
import sys, os, fileinput, string

in_file = open(sys.argv[1], 'r')
in_sam = Reader(in_file)
out_file = open('groundTruth.txt', 'w')
# out_sam = Writer(out_file)

x = next(in_sam)

while(x.qname != ''):
	#if(x.reverse):
	#	out_file.write("+" + " ")
	#else:
	#	out_file.write("-" + " ")
	out_file.write(x.rname + " ")
	out_file.write(x.qname + " ")
	out_file.write(str(x.pos) + " ")
	out_file.write(str(x.pos + len(x.seq)) + "\n")
	#print str(type(x))
	x = next(in_sam)

in_file.close()
out_file.close()