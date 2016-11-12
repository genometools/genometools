#!/bin/sh

seed_extend()
{
  env -i ${GTDIR}/bin/gt seed_extend -display seqlength -minidentity 70 -extend${mode} -ii db -qii query -maxfreq 20 -kmerfile no $1 > sfx-$2.matches
}

# set GTDIR as path of genometools directory

for dbseq in `${GTDIR}/scripts/findfasta.rb`
do
  ${GTDIR}/bin/gt encseq encode -indexname db $dbseq
  for queryseq in `${GTDIR}/scripts/findfasta.rb`
  do
    if test $queryseq != $dbseq
    then
      ${GTDIR}/bin/gt encseq encode -indexname query $queryseq
      for mode in greedy
      do
        echo $filename
        seed_extend "-splt ulong" ulong
        seed_extend "-splt struct" struct
        cmp -s sfx-ulong.matches sfx-struct.matches
        wc -l sfx-ulong.matches
      done
    fi
  done
done
