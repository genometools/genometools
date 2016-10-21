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
    env -i ${GTDIR}/bin/gt -j 4 seed_extend -display seed seqlength -parts 2 -extend${mode} -ii sfx -maxfreq 20 -kmerfile no > sfx-struct.matches
    env -i ${GTDIR}/bin/gt -j 4 seed_extend -display seed seqlength -parts 2 -extend${mode} -ii sfx -maxfreq 20 -kmerfile no -splt ulong > sfx-ulong.matches
    cmp -s sfx-ulong.matches sfx-struct.matches
  done
done
