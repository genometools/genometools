#!/bin/sh

suffixerator()
{
  bin/gt suffixerator -tis -sat ${rep} -db ${1} -indexname esa
}

allmodes()
{
  for mode in stream_words stream_single reader_single stream_reader_single
  do
    bin/gt dev sfxmap -stream-esq esa ${mode}
  done
  for kmersize in 2 3 4 5 6 7 8 9 10 11 12 13
  do
    for mode in reader_multi stream_reader_multi stream_reader_multi3
    do
      bin/gt dev sfxmap -stream-esq esa ${mode} ${kmersize}
    done
  done
}

set -e -x

for filename in `ls testdata/[A-Z]*.fna`
do
  if test ${filename} != "testdata/corruptpatternfile.fna"
  then
    for rep in bit uchar
    do
      suffixerator ${filename}
      allmodes
    done
  fi
done

suffixerator testdata/test1.fastq
allmodes
suffixerator testdata/fastq_long.fastq
allmodes
