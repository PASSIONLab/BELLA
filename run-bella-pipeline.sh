#!/bin/bash
echo "This is a shell script to run BELLA's pipeline and compare accuracy"

# modify as needed
MYPATH=${SCRATCH}/israt-kmer/BELLA

# executable
BELLA=${MYPATH}/./bella
BENCH=${MYPATH}/bench/./result

# input files (modify this as needed)
INPUT=${MYPATH}/input.txt
TRUTH=${MYPATH}/dataset/ecsample-gt.txt

# run parameter (make this input parameter)
DEPTH=30
XDROP=7
KSIZE=17
ERROR=0.15

# todo: implement lower/upper as input parameter
LOWER=2
UPPER=8

# if window is defined and greater than 0, BELLA activates the minimizer counter
WINDOW=0
if [ ${WINDOW} != "0" ]; then
	MMER=true
else
	MMER=false
fi

# syncmer is always false for now, modify the script once syncmer is implemented
SMER=false

if [ ! -d ${MYPATH}/results ]; then
  mkdir -p ${MYPATH}/results;
fi

# choose the name based on the run setting
if   [ $MMER == true ]; then
	NAME="${MYPATH}/results/ecsample-minimizer-${KSIZE}-${WINDOW}-${LOWER}-${UPPER}"
elif [ $SMER == true ]; then	
	NAME="${MYPATH}/results/ecsample-syncmer-${KSIZE}-${WINDOW}-${LOWER}-${UPPER}"
else
	NAME="${MYPATH}/results/ecsample-${KSIZE}-${LOWER}-${UPPER}"
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
echo "	window: ${WINDOW}"
echo "	minimizer: ${MMER}"
echo "	synchmer: ${SMER}"

# creating file for summary result if doesn't exist already
if [ -s ${SUMMARY} ];then
	echo "${SUMMARY} is not empty"
else
	NOW=$(date +"%m-%d-%Y")
	echo $NOW >> ${SUMMARY}
	echo "input	ksize	window	minimizer	syncmer	lower	upper	colA	nnzA	nnzC	nnzR	runtime	recall	precision" >> ${SUMMARY}
fi

MYTEMP="${SCRATCH}/israt-kmer/BELLA/pipeline-tmp-summary.txt"
touch ${MYTEMP}

# run bella (need input for minimizer and synchmer)
${BELLA} -f ${INPUT} -o ${NAME} -c ${DEPTH} -q >> ${MYTEMP}

echo "BELLA run completed"

# collect data from temp file
MYTIME=$(awk 'NR==8 {print; exit}' ${MYTEMP})
TOTKMR=$(awk 'NR==4 {print; exit}' ${MYTEMP})

NNZA=$(awk 'NR==5 {print; exit}' ${MYTEMP})
NNZC=$(awk 'NR==6 {print; exit}' ${MYTEMP})
NNZR=$(awk 'NR==7 {print; exit}' ${MYTEMP})
echo "BELLA data collection completed"

echo "BELLA evaluation is starting"
${BENCH} -G ${TRUTH} -B ${OUTPUT} >> ${MYTEMP}
echo "BELLA evaluation completed"

# todo extract recall and precision
RECALL=
PRECISON=

echo "ecsample	${KSIZE}	${WINDOW}	${MMER}	${SMER}	${LOWER}	${UPPER}	${TOTKMR}	${NNZA}	${NNZC}	${NNZR}	${MYTIME}	${RECALL}	${PRECISION}" >> ${SUMMARY}

# remove tmp summary
# rm ${MYTEMP}

echo "BELLA pipeline completed, results so far can be found here: ${SUMMARY}"

