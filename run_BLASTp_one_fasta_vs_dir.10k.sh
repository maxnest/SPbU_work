#!/bin/bash

QUERY=$1 # fasta file 
QUERY_TAG=$2 # TAG for query fasta file (without .fasta)
DIR_WITH_FASTA=$3 # directory with fasta files
THREADS=$4
EVALUE=$5

### MAIN ###
for SUBJECT in $DIR_WITH_FASTA/*.fasta; do
	SUBJECT_TAG=$(basename $SUBJECT .fasta)
	echo "***** BLASTp starts to work with $QUERY_TAG and $SUBJECT_TAG *****"
	echo "***** e-values = $EVALUE; threads = $THREADS *****"
	echo "***** num_descriptions and num_alignments = 10000 *****"
	### DATABASE CREATING ###
	mkdir ${SUBJECT_TAG}_BLASTdb
	cd ${SUBJECT_TAG}_BLASTdb
	makeblastdb -in $SUBJECT -input_type fasta -dbtype prot -title ${SUBJECT_TAG} -out ${SUBJECT_TAG}_BLASTdb
	wait
	cd ..
	### BLASTp SEARCHING ###
	nohup blastp -db ${SUBJECT_TAG}_BLASTdb/${SUBJECT_TAG}_BLASTdb -query $QUERY -evalue $EVALUE -num_threads $THREADS -num_descriptions 10000 -num_alignments 10000 -out ${QUERY_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}.tab
	wait
	### BBH SEARCHING ###
	echo "***** Best BLAST Hits (BBH) searching between $QUERY_TAG and $SUBJECT_TAG *****"
	nohup /usr/bin/perl /storage/maxnestlcl/Script/parse_best_blast_hit ${QUERY_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}.tab > ${QUERY_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}.bbh.tab
	wait
	nohup rm ${QUERY_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}.tab
	wait
done

