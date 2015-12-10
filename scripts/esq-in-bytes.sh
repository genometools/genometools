#!/bin/sh

set -e -x
# set GTDIR as paht of genometools directory

for filename in `${GTDIR}/scripts/findfasta.rb`
do
  ${GTDIR}/bin/gt encseq encode -indexname sfx $filename
  # Now do something with the sequence
  for mode in greedy xdrop
  do
    ${GTDIR}/bin/gt seed_extend -seedlength 14 -extend${mode} -ii sfx -v -maxfreq 20 -seed-display > sfx.matches
    ${GTDIR}/bin/gt dev show_seedext -f sfx.matches -a
  done
done
