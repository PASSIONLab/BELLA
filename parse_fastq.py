import sys

fastq = sys.argv[1]

with open(fastq) as f:
	content = f.readlines()

	content = [line for line in content if line[0] in 'ACGT']

	data = dict(zip(content[0::2], content[1::2]))

	for x in data:
		print(x) 

