#!/bin/bash
echo "This is a shell script to run BELLA's pipeline and compare accuracy"

# modify as needed
MYPATH=${SCRATCH}/israt-kmer/BELLA

# executable
BELLA=${MYPATH}/./bella
BENCH=${MYPATH}/benchmark/./result

# input files (modify this as needed)
INPUT=${MYPATH}/input.txt
TRUTH=${MYPATH}/dataset/ecsample-truth.txt

# run parameter (make this input parameter)
DEPTH=30 	# dataset-specific
ERROR=0.15 	# dataset-specific

XDROP=$5 # fifth arg

KSIZE=$1 # first arg

# todo: implement lower/upper as input parameter
LOWER=$3 # third arg
UPPER=$4 # fourth arg

# if window is defined and greater than 0, BELLA activates the minimizer counter
WINDOW=$2 # second arg (need to be specific, even if 0)
if [ ${WINDOW} != "0" ]; then
	MMER=true
else
	MMER=false
fi

# syncmer is always false for now, modify the script once syncmer is implemented
SMER=false

if [ ! -d ${MYPATH}/output ]; then
  mkdir -p ${MYPATH}/output;
fi

# choose the name based on the run setting
if   [ $MMER == true ]; then
	NAME="${MYPATH}/output/ecsample-minimizer-${KSIZE}-${WINDOW}-${LOWER}-${UPPER}"
elif [ $SMER == true ]; then	
	NAME="${MYPATH}/output/ecsample-syncmer-${KSIZE}-${WINDOW}-${LOWER}-${UPPER}"
else
	NAME="${MYPATH}/output/ecsample-${KSIZE}-${LOWER}-${UPPER}"
fi

FORMAT=".out"
OUTPUT="${NAME}${FORMAT}"
touch ${OUTPUT}

# one file per dataset
SUMMARY="${MYPATH}/bella-pipeline-ecsample-summary.csv"
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
	echo "input	xdrop	ksize	window	minimizer	syncmer	lower	upper	colA	nnzA	nnzC	nnzR	time	RC	PR	F1" >> ${SUMMARY}
fi

MYTEMP="${SCRATCH}/israt-kmer/BELLA/pipeline-tmp-summary.txt"
touch ${MYTEMP}

if [ ${WINDOW} == "0" ]; then
	# run bella (need input for minimizer and synchmer)
	${BELLA} -f ${INPUT} -k ${KSIZE} -o ${NAME} -x ${XDROP} -l ${LOWER} -u ${UPPER} -t >> ${MYTEMP}
else
	${BELLA} -f ${INPUT} -k ${KSIZE} -o ${NAME} -w ${WINDOW} -x ${XDROP} -l ${LOWER} -u ${UPPER} -t >> ${MYTEMP}
fi

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
RC=$(awk 'NR==16 {print; exit}' ${MYTEMP})
PR=$(awk 'NR==17 {print; exit}' ${MYTEMP})
F1=$(awk 'NR==18 {print; exit}' ${MYTEMP})

echo "ecsample	${XDROP}	${KSIZE}	${WINDOW}	${MMER}	${SMER}	${LOWER}	${UPPER}	${TOTKMR}	${NNZA}	${NNZC}	${NNZR}	${MYTIME}	${RC}	${PR}	${F1}" >> ${SUMMARY}

# remove tmp summary
rm ${MYTEMP}

echo "BELLA pipeline completed, output so far can be found here: ${SUMMARY}"

