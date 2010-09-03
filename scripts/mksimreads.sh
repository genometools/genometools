#!/bin/sh

if test $# -ne 3
then
  echo "Usage: $0 <readsize> <coverage> <inputfile>"
  exit 1
fi

readsize=$1
coverage=$2
inputfile=$3

gt suffixerator -des -tis -ssp -dna -db ${inputfile} -indexname genome-idx
gt simreads -coverage ${coverage} -len ${readsize} -gzip -force -o genome-idx-reads.fna.gz genome-idx
rm -f genome-idx.*
