#!/bin/sh

set -e -x

for filename in `ls testdata/[A-Z]*.fna`
do
  for rep in bit uchar
  do
    bin/gt suffixerator -tis -sat ${rep} -db ${filename} -indexname esa
    bin/gt dev sfxmap -stream-esq esa
  done
done

bin/gt suffixerator -tis -db testdata/test1.fastq -indexname esa
bin/gt dev sfxmap -stream-esq esa
bin/gt suffixerator -tis -db testdata/fastq_long.fastq -indexname esa
bin/gt dev sfxmap -stream-esq esa
