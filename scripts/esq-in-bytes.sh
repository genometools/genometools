#!/bin/sh

set -e -x
# set GTDIR as paht of genometools directory

for filename in `${GTDIR}/scripts/findfasta.rb`
do
  len=`cat ${filename} | wc -c`
  if test ${len} -ge 20000
  then
    ${GTDIR}/bin/gt encseq encode -indexname sfx $filename
    # Now do something with the sequence
    ${GTDIR}/bin/gt seed_extend -extendgreedy -ii sfx -a -maxfreq 20
  fi
done
