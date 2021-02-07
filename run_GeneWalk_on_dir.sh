#!/bin/bash

INPUT_DIR=$1
THREADS=$2

### ENV ###

GeneWalk_env=/sf/smpdata1/maxnestlcl/GeneWalk_env/genewalk_env/bin

### MAIN ###
cd $INPUT_DIR

source $GeneWalk_env/activate

for GeneSet in $INPUT_DIR/*.txt; do
	GeneSet_TAG=$(basename $GeneSet .txt)
	echo "***** The GeneWalk starts to work with $GeneSet_TAG *****"
	### OUTPUT DIR CREATING ###
	mkdir ${INPUT_DIR}/${GeneSet_TAG}_genewalk_output
	wait
	### ANALYSIS ###
	nohup genewalk --project $GeneSet_TAG --base_folder ${INPUT_DIR}/${GeneSet_TAG}_genewalk_output --genes $GeneSet --id_type ensembl_id --nproc $THREADS
	wait
	echo "***** The GeneWalk complete analysis for $GeneSet_TAG *****"
done

deactivate
