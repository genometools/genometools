#!/bin/sh
#
# Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

# Current Problems: for icc compiled version add 
# /usr/local/zbh/intel/cc/9.1/lib

set -e -x

outoptions="-tis -lcp -suf -bwt"

# the make call normally used for development
cd testsuite
testsuite.rb -keywords gt_suffixerator
# optional -memcheck   (run valgrind)
#          -select 253 (run testcase 253)
../scripts/runmaxpairs.sh 14 ${GRUMBACH}/*.fna
num=2
while test ${num} -lt 10 
do
  ../scripts/iterrunmerge.sh ${num}
  num=`expr ${num} + 1`
done
for options in `alloutputoptions.rb`
do
  cmpdbfile.sh ${options} -pl -db ${ATK}
done
../scripts/checkmapped.sh -db ${ATK} ${AT} ${GRUMBACH}/*.fna -parts 3 -pl
../scripts/checkmapped.sh -parts 1 -pl -db ${SWK} ${SW}
../scripts/checkmapped.sh -db ${SWK} ${SW} -parts 3 -pl
../scripts/checkmapped.sh -parts 2 -pl -smap TransDNA -db ${AT}
../scripts/checkmapped.sh -db ${SWK} -parts 1 -pl -smap TransProt11 
../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/Random-Small.fna
../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/Random.fna
../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/Atinsert.fna ../testdata/Random.fna
../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/TTT-small.fna
../scripts/cmpdbfile.sh ${outoptions} -pl -db ${ATK} ${AT} ${GRUMBACH}/*.fna
testsuite.rb -keywords gt_trieins
cd ..
