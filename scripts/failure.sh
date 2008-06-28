#!/bin/sh

set -e -x

# exec in gt root dir
# works without -sprank and fails with -sprank

#SPRANK="-sprank"

bin/gt packedindex mkindex ${SPRANK} -tis -dna -pl -bsize 10 -locfreq 32\
                                     -dir rev -db ${U8}

# bin/gt tagerator -t Q1 -pck pck
