#!/bin/sh
bsize=4

gt packedindex mkindex -protein -pl -tis -ssp -sprank -bsize ${bsize} \
                       -locfreq 128 -dir rev -parts 20 -v \
                       -db /local/kurtz/sfx-test/data/trembl-fasta.gz \
                       -indexname tremb64-${bsize}
