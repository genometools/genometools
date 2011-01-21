#!/bin/sh

suffixerator()
{
  bin/gt suffixerator -tis -sat ${rep} -db ${1} -indexname esa
}

allmodes()
{
  for mode in 0 1 2
  do
    bin/gt dev sfxmap -stream-esq esa ${mode}
  done
}

set -e -x

for filename in `ls testdata/[A-Z]*.fna`
do
  for rep in bit uchar
  do
    suffixerator ${filename}
    allmodes
  done
done

suffixerator testdata/test1.fastq
allmodes
suffixerator testdata/fastq_long.fastq
allmodes
