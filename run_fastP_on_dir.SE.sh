#!/bin/bash

# INPUT : #
data_dir=$1

# SOFT: #

FASTP=/storage/maxnestlcl/SOFT/fastp/fastp

# MAIN: #
for dir in $(find $data_dir -mindepth 1 -type d); do
	lib="$(find $dir -type f -name '*.fastq.gz')"
        tag=$(basename $dir)
        echo "The fastP began to work with params: "
        echo "Lib: $lib"
        echo "TAG: $tag"
        mkdir ${tag}_fastp_out
        cd ${tag}_fastp_out
	nohup $FASTP -i $lib -o ${tag}.fastp_tmm.fastq.gz --failed_out ${tag}.fastp_failed_out.fastq.gz --trim_poly_x --cut_right --cut_window_size 4 --cut_mean_quality 20 --qualified_quality_phred 20 --length_required 25 --html ${tag}.fastp_report.html
	echo "### ${tag} complete ###"
	cd ..
done
