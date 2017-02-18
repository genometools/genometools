#!/bin/sh

# Take <numreads> samples from <inputfile> and map the samples to the reference. 
# Use collect-mappings.rb for determining sensitivity and save it to a table.

set -e

if test $# -ne 2
then
  echo "Usage: $0 <numreads> <inputfile>"
  exit 1
fi
numreads=$1
inputfile=$2
readlength=150
readset=reads.fa
outfile_nw_m=greedy-noweak-mye.tsv
outfile_w_m=greedy-weak-mye.tsv
outfile_nw_nm=greedy-noweak-nomye.tsv
outfile_w_nm=greedy-weak-nomye.tsv
gt_bin="./bin"
#da_bin="../myers/"

minlength=`expr ${readlength} \* 80`
minlength=`expr ${minlength} \/ 100`
$gt_bin/gt encseq encode -indexname reference -sds no -md5 no -des no ${inputfile}
#$gt_bin/../scripts/convert2myersformat.rb $inputfile > reference.fasta
#$da_bin/DAZZ_DB/fasta2DB reference.db reference.fasta
#inputfile=reference.fasta

for minid in {70..99}; do printf "\t%d" $minid; done > $outfile_nw_m
for minid in {70..99}; do printf "\t%d" $minid; done > $outfile_w_m
for minid in {70..99}; do printf "\t%d" $minid; done > $outfile_nw_nm
for minid in {70..99}; do printf "\t%d" $minid; done > $outfile_w_nm

for err in {0..30}
do
  printf "\n%d" $err >> $outfile_nw_m
  printf "\n%d" $err >> $outfile_w_m
  printf "\n%d" $err >> $outfile_nw_nm
  printf "\n%d" $err >> $outfile_w_nm

  indel=$(echo "scale=6; ${err}/3000"|bc)
  misma=$(echo "scale=6; ${err}/100-2*${indel}"|bc)
  mason_simulator -ir $inputfile -n $numreads -o ${readset} \
                --illumina-read-length ${readlength} --embed-read-info \
                --illumina-prob-mismatch $misma \
                --illumina-prob-deletion $indel \
                --illumina-prob-insert $indel
  $gt_bin/gt encseq encode -indexname query -sds no -md5 no -des no ${readset}
  #$gt_bin/../scripts/convert2myersformat.rb $readset > query.fasta
  #$da_bin/DAZZ_DB/fasta2DB query.db query.fasta

  for minid in {70..99}
  do
    common="-l ${minlength} -v -seedlength 14 -ii reference -qii query -minidentity $minid"
    $gt_bin/gt seed_extend ${common} -bias-parameters > tmp.matches
    #$da_bin/DALIGNER/daligner -l${minlength} -I -A -Y -k14 -t21 -e0.${minid} \
    #                          query.db reference.db > tmp.matches
    #grep '^#' tmp.matches
    sensitivity=`./collect-mappings.rb ${numreads} tmp.matches | grep -v '^#'`
    printf "\t${sensitivity}" >> $outfile_nw_m
    rm -f tmp.matches *.las

    $gt_bin/gt seed_extend ${common} -weakends -bias-parameters > tmp.matches
    sensitivity=`./collect-mappings.rb ${numreads} tmp.matches | grep -v '^#'`
    printf "\t${sensitivity}" >> $outfile_w_m
    rm -f tmp.matches *.las

    $gt_bin/gt seed_extend ${common} > tmp.matches
    sensitivity=`./collect-mappings.rb ${numreads} tmp.matches | grep -v '^#'`
    printf "\t${sensitivity}" >> $outfile_nw_nm
    rm -f tmp.matches *.las

    $gt_bin/gt seed_extend ${common} -weakends > tmp.matches
    sensitivity=`./collect-mappings.rb ${numreads} tmp.matches | grep -v '^#'`
    printf "\t${sensitivity}" >> $outfile_w_nm
    rm -f tmp.matches *.las
  done
  rm -f ${readset} query.*
done
rm -f reference.*
