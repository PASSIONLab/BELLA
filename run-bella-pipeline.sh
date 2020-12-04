#!/bin/bash
echo "This is a shell script to run BELLA's pipeline and compare accuracy"

# executable (change this as needed)
BELLA=$SCRATCH/israt-kmer/BELLA/./bella
BENCH=$SCRATCH/israt-kmer/BELLA/bench/./result

# input and output files
INPUT=
OUTPUT=
SUMMARY=
TRUTH=

# run parameter (make this input parameter)
DEPTH=30
XDROP=7
KSIZE=17
ERROR=0.15
MINIMIZER=true
SYNCMER=false

# start run bella
echo "${now}"
echo "BELLA run is starting:"
echo "	input: 	${INPUT}"
echo "	output: ${OUTPUT}"
echo "	depth: 	${DEPTH}"
echo "	x-drop: ${XDROP}"
echo "	k-mer:  ${KSIZE}"
echo "	error:  ${ERROR}"
echo "  minimizer: ${MINIMIZER}"
echo "	synchmer: ${SYNCMER}"

# creating file for summary result if doesn't exist already
if [-s ${SUMMARY}]
then
	echo "${SUMMARY} is not empty"
else
	echo "${now}"
	echo "input\tksize\twindow\tminimizer\tsyncmer\truntime\trecall\tprecision\tnalignment" >> ${SUMMARY}
fi

MYTEMP=${SCRATCH}/israt-kmer/BELLA/pipeline-tmp-summary.txt

# run bella (need input for minimizer and synchmer)
${BELLA} -f ${INPUT} -o ${OUTPUT} -c ${DEPTH} -q >> ${MYTEMP}

echo "BELLA run completed"

# GGGG: modify main.cpp to make the time retrival easy
# todo extract runtime and nalignment
MYTIME=
NALIGN=

echo "BELLA evaluation is starting"

${BENCH} -G ${TRUTH} -B ${OUTPUT} >> ${MYTEMP}

echo "BELLA evaluation completed"

# todo extract recall and precision
RECALL=
PRECIISON=

echo "${INPUT}\t${KSIZE}\t${WINDOW}\t${MINIMZER}\t${SYNCMER}\t${MYTIME}\t${RECALL}\t${PRECISION}" >> ${SUMMARY}

# remove tmp summary
rm ${MYTEMP}

echo "BELLA pipeline completed, results so far can be found here: ${SUMMARY}"

