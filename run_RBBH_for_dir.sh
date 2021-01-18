#!/bin/bash

DIR_WITH_BBH=$1 

### SCRIPTS ###
RBBH_SCRIPT=/storage/maxnestlcl/Script/RBBH_pairs.py

### MAIN ###
for file in $DIR_WITH_BBH/*.bbh.tab; do
	first_sp=$(basename $file | awk -F "_vs_" '{print $1}')
	second_sp=$(basename $file | awk -F "_vs_" '{print $2}' | awk -F "_10k_" '{print $1}')
	second_file=${second_sp}_vs_${first_sp}_10k_1e-3.bbh.tab
	cd ${DIR_WITH_BBH}
	# echo "***** RBBH searching between $first_sp and $second_sp *****"
	python3.6 $RBBH_SCRIPT --blast_sp1_vs_sp2 $file --blast_sp2_vs_sp1 ${DIR_WITH_BBH}/$second_file --sp1_tag $first_sp --sp2_tag $second_sp
	wait
done
