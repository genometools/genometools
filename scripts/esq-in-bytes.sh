#!/bin/sh

set -e -x
# set GTDIR as paht of genometools directory

for filename in `${GTDIR}/scripts/findfasta.rb`
do
  gt encseq encode -ssp no -sds no -md5 no -des no -no_esq_header \
                   -showstats -sat direct -indexname seq $filename
  # Now do something with the sequence
  ls -l seq.esq
done
