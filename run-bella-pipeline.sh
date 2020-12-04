#!/bin/bash
echo "This is a shell script to run BELLA's pipeline and compare accuracy"

# executable (change this as needed)
BELLA=${SCRATCH}/israt-kmer/BELLA/./bella
BENCH=${SCRATCH}/israt-kmer/BELLA/bench/./result

# input files (change this as needed)
INPUT=${SCRATCH}/israt-kmer/BELLA/input.txt
TRUTH=${SCRATCH}/israt-kmer/BELLA/dataset/ecsample-gt.txt

# run parameter (make this input parameter)
DEPTH=30
XDROP=7
KSIZE=17
ERROR=0.15

# todo: implement lower/upper as input parameter
LOWER=
UPPER=

# if window is defined and greater than 0, BELLA activates the minimizer counter
WINDOW=0
if [ $WINDOW > 0]; then
	MMER=true
else
	MMER=false

# syncmer is always false for now, modify the script once syncmer is implemented
SMER=false

# choose the name based on the run setting
if   [ $MMER == true ]; then
	NAME="${SCRATCH}/israt-kmer/BELLA/ecsample-minimizer-${KSIZE}-${WINDOW}"
elif [ $SMER == true ]; then	
	NAME="${SCRATCH}/israt-kmer/BELLA/ecsample-syncmer-${KSIZE}-${WINDOW}"
else
	NAME="${SCRATCH}/israt-kmer/BELLA/ecsample-${KSIZE}"
fi

FORMAT=".out"
OUTPUT="${NAME}${FORMAT}"
touch ${OUTPUT}

SUMMARY="${SCRATCH}/israt-kmer/BELLA/bella-pipeline-summary.csv"
touch ${SUMMARY}

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
if [ -s ${SUMMARY} ]
then
	echo "${SUMMARY} is not empty"
else
	NOW=$(date +"%m-%d-%Y")
	echo $NOW >> ${SUMMARY}
	echo "input	ksize	window	minimizer	syncmer	runtime	recall	precision	nalignment" >> ${SUMMARY}
fi

MYTEMP=${SCRATCH}/israt-kmer/BELLA/pipeline-tmp-summary.txt
touch ${MYTEMP}

# run bella (need input for minimizer and synchmer)
${BELLA} -f ${INPUT} -o ${NAME} -c ${DEPTH} -q >> ${MYTEMP}

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

echo "ecsample	${KSIZE}	${WINDOW}	${MINIMIZER}	${SYNCMER}	${MYTIME}	${RECALL}	${PRECISION}" >> ${SUMMARY}

# remove tmp summary
$rm ${MYTEMP}

echo "BELLA pipeline completed, results so far can be found here: ${SUMMARY}"

