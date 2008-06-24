#!/bin/sh

set -e -x

# exec in gt root dir
# works without -sprank and fails with -sprank

SPRANK="-sprank"

bin/gt packedindex mkindex ${SPRANK} -tis -dna -pl -bsize 10 -locfreq 32\
                                     -dir rev -indexname pck \
                                     -db testdata/Atinsert.fna

bin/gt matstat -verify -output querypos -min 1 -max 20\
               -query testdata/Duplicate.fna -pck pck
