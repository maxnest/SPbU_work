#!/bin/bash

focal_fasta=$1
uniprot_seqs=$2

### MAIN ###

for uniprot_fasta in $(find $uniprot_seqs -type f); do
    tag=$(basename $uniprot_fasta)
    echo "Comparison with $uniprot_fasta"
    mkdir ${tag}_diamond_blast_out
    cd ${tag}_diamond_blast_out
    nohup diamond makedb --in $uniprot_fasta -d $tag
    wait
    echo -e 'qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tevalue\tscore' > ${tag}.tab
    nohup diamond blastp -q $focal_fasta -d $tag --more-sensitive --threads 5 --outfmt 6 qseqid sseqid qstart qend sstart send evalue score >> ${tag}.tab
    wait
    cd ..
done
