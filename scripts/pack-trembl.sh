#!/bin/sh

set -e -x

# INPUTFILE=/local/kurtz/sfx-test/data/trembl-fasta.gz
INPUTFILE=${SWISSDIR}/swissprot
BINDIR=./bin
WSIZE=64
locfreq=64

for bsize in 1
do
  ${BINDIR}/gt packedindex mkindex -protein -pl -tis -ssp -sprank \
                                   -bsize ${bsize} \
                                   -locfreq ${locfreq} -dir rev -parts 20 -v \
                                   -db ${INPUTFILE} \
                                   -indexname pck${WSIZE}-${bsize} >\
                                   ${HOME}/PROTO${WSIZE}-${bsize}.txt
done
