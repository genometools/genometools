#!/bin/sh
#
# Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

# Current Problems: for icc compiled version add 
# /usr/local/zbh/intel/cc/9.1/lib

set -e -x

outoptions="-tis -lcp -suf -bwt"
ALLOUTPUTOPTS="../scripts/alloutputoptions.rb"

# the make call normally used for development
cd testsuite
testsuite.rb -keywords gt_suffixerator
testsuite.rb -keywords gt_trieins
# optional -memcheck   (run valgrind)
#          -select 253 (run testcase 253)
# the following depends on vmatch-mini.x and mkvtree.x
# ../scripts/runmaxpairs.sh 14 ${GRUMBACH}/*.fna ../testdata/Duplicate.fna
num=2
while test ${num} -lt 10 
do
  ../scripts/iterrunmerge.sh ${num}
  num=`expr ${num} + 1`
done
if test ! -f ${ALLOUTPUTOPTS}
then
  echo "cannot find ${ALLOUTPUTOPTS}"
  exit 1
fi
for options in `${ALLOUTPUTOPTS}`
do
  ../scripts/cmpdbfile.sh ${options} -pl -db ${ATK}
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
cd ..
