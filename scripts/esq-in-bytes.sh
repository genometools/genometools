#!/bin/sh

set -e -x
# set GTDIR as paht of genometools directory

for filename in `${GTDIR}/scripts/findfasta.rb`
do
  len=`cat ${filename} | wc -c`
  ${GTDIR}/bin/gt encseq encode -indexname sfx $filename
  # Now do something with the sequence
  ${GTDIR}/bin/gt seed_extend -extendgreedy -ii sfx -a -maxfreq 20 -seed-display
done
