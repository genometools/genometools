#!/bin/sh

set -e -x

DATA="/local/kurtz/sfx-test/data"

run()
{
  gt suffixerator $1 -tis -suf -pl -showtime -db $2
}

for filename in `ls ${DATA}/fib*.fas.gz`
do
  run "-smap Transab" "${filename}"
done
#run "-protein" "${DATA}/uniprot_sprot.fas.gz"
#run "-dna" "`humanpath.sh 1`"
#run "-dna" "${HOME}/seqcmpprojects/MouthFootDisease/mfdall.fna.gz"
