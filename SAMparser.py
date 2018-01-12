from simplesam import Reader, Writer, Sam
import inspect
import sys, os, fileinput, string

in_file = open(sys.argv[1], 'r')
in_sam = Reader(in_file)
out_file = open('PBcR-PB-ec-Q10-strand.txt', 'w')
# out_sam = Writer(out_file)

x = next(in_sam)

while(x != ''):
	if(x.reverse):
		out_file.write("+" + " ")
	else:
		out_file.write("-" + " ")
	out_file.write(x.qname + " ")
	out_file.write(str(x.pos) + " ")
	out_file.write(str(x.pos + len(x.seq)) + "\n")
	#print str(type(x))
	x = next(in_sam)

in_file.close()
out_file.close()
