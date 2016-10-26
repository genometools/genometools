#!/bin/sh

seed_extend()
{
  env -i ${GTDIR}/bin/gt seed_extend -display seed seqlength -parts 2 -extend${mode} -ii sfx -maxfreq 20 -kmerfile no $1 > sfx-$2.matches
}

set -e -x
# set GTDIR as path of genometools directory

for filename in `${GTDIR}/scripts/findfasta.rb`
do
  ${GTDIR}/bin/gt encseq encode -indexname sfx $filename
  # Now do something with the sequence
  for mode in greedy xdrop
  do
    echo $filename
    seed_extend "-splt ulong" ulong
    cmp -s sfx-ulong.matches sfx-struct.matches
    seed_extend "-splt struct" bytestring
    cmp -s sfx-ulong.matches sfx-bytestring.matches
  done
done
