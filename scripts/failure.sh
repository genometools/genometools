#!/bin/sh

set -e -x

if test $# -ne 2
then
  echo "Usage: $0 <Referencefile> <Queryfile>"
  exit 1
fi  

reference=$1
query=$2

# exec in gt root dir
# works without -sprank and fails with -sprank

#SPRANK="-sprank"

gt suffixerator -tis -suf -dna -pl -db ${reference} -indexname esa-fwd
gt packedindex mkindex -sprank -tis -dna -pl -bsize 10 -locfreq 32\
                       -dir rev -db ${reference} -indexname pck-rev
gt tagerator -t ${query} -cmp -esa esa-fwd -k 1
gt tagerator -t ${query} -cmp -pck pck-rev -k 1

# bin/gt tagerator -t Q1 -pck pck
