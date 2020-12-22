QUERY=$1 # fasta file
SUBJECT=$2 # fasta file
INTER=$3 # fasta file
THREADS=$4
EVALUE=$5
QUERY_TAG=$6 
SUBJECT_TAG=$7
INTER_TAG=$8

##### MAIN #####

echo "***** STEP 0:  DATABASES creating *****"
### New directories ###
mkdir ${QUERY_TAG}_BLASTDB
mkdir ${SUBJECT_TAG}_BLASTDB
mkdir ${INTER_TAG}_BLASTDB
### makeblastdb ###
## 1 ##
cd ${QUERY_TAG}_BLASTDB
makeblastdb -in $QUERY -input_type fasta -dbtype prot -title ${QUERY_TAG} -out ${QUERY_TAG}_BLASTDB
wait
cd ..
## 2 ##
cd ${SUBJECT_TAG}_BLASTDB
makeblastdb -in $SUBJECT -input_type fasta -dbtype prot -title ${SUBJECT_TAG} -out ${SUBJECT_TAG}_BLASTDB
wait
cd ..
## 3 ## 
cd ${INTER_TAG}_BLASTDB
makeblastdb -in $INTER -input_type fasta -dbtype prot -title ${INTER_TAG} -out ${INTER_TAG}_BLASTDB
wait
cd ..

### BLAST searching ###
## 1 vs 2 ##
echo "***** Step 1: BLASTp starts to work with $QUERY_TAG and $INTER_TAG (e-values = $EVALUE; outfmt = 6) *****"
nohup blastp -query $QUERY -db ${INTER_TAG}_BLASTDB/${INTER_TAG}_BLASTDB -evalue $EVALUE -num_threads $THREADS -max_target_seqs 10000 -outfmt 6 -out ${QUERY_TAG}_vs_${INTER_TAG}_10k_${EVALUE}.tab
wait
echo "***** STEP 2: Best BLAST Hits (BBH) searching between $QUERY_TAG and $INTER_TAG *****"
nohup python3.6 /storage/maxnestlcl/Script/DIAMOND_outfmt6_bbh.py --outfmt_6 ${QUERY_TAG}_vs_${INTER_TAG}_10k_${EVALUE}.tab --output ${QUERY_TAG}_vs_${INTER_TAG}_10k_${EVALUE}
wait
nohup rm ${QUERY_TAG}_vs_${INTER_TAG}_10k_${EVALUE}.tab
wait

## 2 vs 1 ##
echo "***** STEP 3: BLASTp starts to work with $INTER_TAG and $QUERY_TAB (e-value = $EVALUE; outfmt = 6) *****"
nohup blastp -query $INTER -db ${QUERY_TAG}_BLASTDB/${QUERY_TAG}_BLASTDB -evalue $EVALUE -num_threads $THREADS -max_target_seqs 10000 -outfmt 6 -out ${INTER_TAG}_vs_${QUERY_TAG}_10k_${EVALUE}.tab
wait
echo "***** STEP 4: Best BLAST Hits (BBH) searching between $INTER_TAG and $QUERY_TAG *****"
nohup python3.6 /storage/maxnestlcl/Script/DIAMOND_outfmt6_bbh.py --outfmt_6 ${INTER_TAG}_vs_${QUERY_TAG}_10k_${EVALUE}.tab --output ${INTER_TAG}_vs_${QUERY_TAG}_10k_${EVALUE}
wait
nohup rm ${INTER_TAG}_vs_${QUERY_TAG}_10k_${EVALUE}.tab
wait

## 1 vs 3 ##
echo "***** STEP 5: BLASTP starts to work with $QUERY_TAG and $SUBJECT_TAG (e-values = $EVALUE; outfmt = 6) *****"
nohup blastp -query $QUERY -db ${SUBJECT_TAG}_BLASTDB/${SUBJECT_TAG}_BLASTDB -evalue $EVALUE -num_threads $THREADS -max_target_seqs 10000 -outfmt 6 -out ${QUERY_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}.tab
wait
echo "***** STEP 6: Best BLAST Hits (BBH) searching between $QUERY_TAG and $SUBJECT_TAG *****"
nohup python3.6 /storage/maxnestlcl/Script/DIAMOND_outfmt6_bbh.py --outfmt_6 ${QUERY_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}.tab --output ${QUERY_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}
wait
nohup rm ${QUERY_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}.tab
wait

## 2 vs 3 ##
echo "***** STEP 7: BLASTP starts to work with $INTER_TAG and $SUBJECT_TAG (e-values = $EVALUE; outfmt = 6) *****"
nohup blastp -query $INTER -db ${SUBJECT_TAG}_BLASTDB/${SUBJECT_TAG}_BLASTDB -evalue $EVALUE -num_threads $THREADS -max_target_seqs 10000 -outfmt 6 -out ${INTER_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}.tab
wait
echo "***** STEP 8: Best BLAST Hits (BBH) searching between $INTER_TAG and $SUBJECT_TAG *****"
nohup python3.6 /storage/maxnestlcl/Script/DIAMOND_outfmt6_bbh.py --outfmt_6 ${INTER_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}.tab --output ${INTER_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}
wait
nohup rm ${INTER_TAG}_vs_${SUBJECT_TAG}_10k_${EVALUE}.tab
wait

## 3 vs 2 ##
echo "***** STEP 9: BLASTP starts to work with $SUBJECT_TAG and $INTER_TAG (e-values = $EVALUE; outfmt = 6) *****"
nohup blastp -query $SUBJECT -db ${INTER_TAG}_BLASTDB/${INTER_TAG}_BLASTDB -evalue $EVALUE -num_threads $THREADS -max_target_seqs 10000 -outfmt 6 -out ${SUBJECT_TAG}_vs_${INTER_TAG}_10k_${EVALUE}.tab
wait
echo "***** STEP 10: Best BLAST Hits (BBH) searching between $SUBJECT_TAG and $INTER_TAG *****"
nohup python3.6 /storage/maxnestlcl/Script/DIAMOND_outfmt6_bbh.py --outfmt_6 ${SUBJECT_TAG}_vs_${INTER_TAG}_10k_${EVALUE}.tab --output ${SUBJECT_TAG}_vs_${INTER_TAG}_10k_${EVALUE}
wait
nohup rm ${SUBJECT_TAG}_vs_${INTER_TAG}_10k_${EVALUE}.tab
wait
#####
echo "***** Finish *****"
