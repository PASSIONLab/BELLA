#!/bin/bash
echo "This is a shell script to run BELLA's pipeline in batches"

./run-bella-pipeline.sh 13 21 2 10 25
./run-bella-pipeline.sh 14 21 2 10 25
./run-bella-pipeline.sh 15 21 2 10 25

./run-bella-pipeline.sh 13 19 2 10 25
./run-bella-pipeline.sh 14 19 2 10 25
./run-bella-pipeline.sh 15 19 2 10 25

./run-bella-pipeline.sh 17 0 2 7 25

echo "Batch completed"
