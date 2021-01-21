import sys
from Bio import SeqIO
from Bio.Seq import Seq

MAX_READ_LENGTH = 20272

def ErrorCal(fastqfile):

    f = open(fastqfile)

    prob_sums = [0 for i in range(MAX_READ_LENGTH)]
    prob_pos = prob_sums[:]
    read_count = 0
    for record in SeqIO.parse(f, "fastq"):
        for bp, Q in enumerate(record.letter_annotations["phred_quality"]):
            P = 10**(float(-Q)/10)
            prob_sums[bp] += P
            prob_pos[bp] += 1.0

    for bp, data in enumerate(zip(prob_sums, prob_pos)):
        prob_sum, prob_pos = data
        prob_mean = prob_sum / prob_pos
        print(bp, prob_pos, prob_mean)

if __name__ == '__main__':
    ErrorCal(sys.argv[1])