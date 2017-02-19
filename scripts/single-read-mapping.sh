#!/bin/sh

set -e -x

if test $# -ne 4
then
  echo "Usage: $0 <numreads> <inputfile> <mincoverage> <seedlength>"
  exit 1
fi
numreads=$1
inputfile=$2
mincoverage=$3
seedlength=$4
readlength=150
maxfreq=30
minidentity=95
readset=reads.fa
#probabilities='--illumina-prob-mismatch 0.028 --illumina-prob-deletion 0.001 --illumina-prob-insert 0.001'
mason_simulator -ir $inputfile -n $numreads -o ${readset} \
                --illumina-read-length ${readlength} --embed-read-info \
                ${probabilities}

env -i bin/gt encseq encode -indexname query-idx ${readset}
env -i bin/gt encseq encode -indexname reference-idx ${inputfile}
minlength=`expr ${readlength} \* ${mincoverage}`
minlength=`expr ${minlength} \/ 100`
common="-l ${minlength} -v -ii reference-idx -qii query-idx -minidentity ${minidentity} -seedlength ${seedlength}"
env -i bin/gt seed_extend ${common} > tmp.matches
scripts/collect-mappings.rb ${readset} tmp.matches
env -i bin/gt seed_extend ${common} -maxmat 2 -use-apos > tmp-maxmat.matches
scripts/collect-mappings.rb ${readset} tmp-maxmat.matches
