#!/bin/bash

# INPUT : #
data_dir=$1

# SOFT: #

FASTP=/storage/maxnestlcl/SOFT/fastp/fastp

# MAIN: #
for dir in $(find $data_dir -mindepth 1 -type d); do
	r1="$(find $dir -type f -name '*.R1.fastq')"
        r2="$(find $dir -type f -name '*.R2.fastq')"
        tag=$(basename $dir)
        echo "The fastP began to work with params: "
        echo "R1: $r1"
        echo "R2: $r2"
        echo "TAG: $tag"
        mkdir ${tag}_fastp_out
        cd ${tag}_fastp_out
	nohup $FASTP --in1 $r1 --out1 ${tag}.fastp_tmm.R1.fastq --in2 $r2 --out2 ${tag}.fastp_tmm.R2.fastq --unpaired1 ${tag}.fastp_unpaired1.fastq --unpaired2 ${tag}.fastp_unpaired2.fastq --failed_out ${tag}.fastp_failed_out.fastq --trim_poly_x --cut_right --cut_window_size 4 --cut_mean_quality 20 --qualified_quality_phred 20 --length_required 25 --html ${tag}.fastp_report.html
	echo "### ${tag} complete ###"
	cd ..
done
