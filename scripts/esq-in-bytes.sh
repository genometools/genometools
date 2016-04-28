#!/bin/sh

set -e -x
# set GTDIR as path of genometools directory

for filename in `${GTDIR}/scripts/findfasta.rb`
do
  ${GTDIR}/bin/gt encseq encode -indexname sfx $filename
  # Now do something with the sequence
  for mode in greedy xdrop
  do
    echo $filename
    ${GTDIR}/bin/gt -j 4 seed_extend -seqlength-display -parts 2 -extend${mode} -ii sfx -v -maxfreq 20 -seed-display -kmerfile no > sfx.matches
    ${GTDIR}/bin/gt dev show_seedext -seqlength-display -f sfx.matches -a -sort
  done
done
