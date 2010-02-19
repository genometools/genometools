#!/bin/sh
bsize=3

gt suffixerator -protein -pl -tis -ssp -suf -des -lcp -bwt \
                -parts 20 -v \
                -db /local/kurtz/sfx-test/data/trembl-fasta.gz \
                -indexname sarr-tremb32-${bsize}
