#!/bin/sh
Yeast=`ls ${GTTESTDATA}/ltrharvest/s_cer/chr[01][0-9].*.gz`
bin/gt packedindex mkindex -sprank -tis -dna -pl -bsize 10 -locfreq 32 -dir rev\
                       -db ${Yeast}\
                       -indexname pck-yeast

env -i GT_MEM_BOOKKEEPING=off time bin/gt tagerator -rw -pck pck-yeast -k 2\
    -t yeast.1000 > tmp.prot
gprof bin/gt gmon.out
