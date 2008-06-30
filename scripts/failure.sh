#!/bin/sh

set -e -x

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi  

inputfile=$1

# exec in gt root dir
# works without -sprank and fails with -sprank

#SPRANK="-sprank"

gt suffixerator -tis -suf -dna -pl -db ${inputfile} -indexname esa-fwd
gt suffixerator -tis -bwt -suf -dna -pl -db ${inputfile} -dir rev -indexname esa-rev
gt tagerator -t Q1 -cmp -esa esa-fwd -k 1
gt packedindex mkindex -tis -dna -pl -bsize 10 -locfreq 32\
                       -dir rev -db ${inputfile} -indexname pck-rev
gt tagerator -t Q1 -cmp -pck pck-rev -k 1

# bin/gt tagerator -t Q1 -pck pck
