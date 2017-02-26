#!/bin/sh

set -e -x

seed_extend()
{
  env -i ${GTDIR}/bin/gt seed_extend -l 40 -outfmt alignment=50  -minidentity 70  -ii db -qii query -maxfreq 20 -kmerfile no $1
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
        seed_extend "-no-reverse -maxmat 2 -use-apos"
      done
    fi
  done
done
