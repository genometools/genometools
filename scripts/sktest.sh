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

USAGE="Usage: $0 [-memcheck]"

if test $# -eq 0
then
  MC=""
else
  if test $# -eq 1
  then
    if test "$1" = "-memcheck"
    then
      MC="-memcheck" 
    else
      echo ${USAGE}
      exit 1
    fi 
  else
    echo ${USAGE}
    exit 1
  fi
fi

outoptions="-tis -lcp -suf -bwt"
ALLOUTPUTOPTS="../scripts/alloutputoptions.rb"

cd testsuite

# the make call normally used for development
env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb ${MC} -keywords 'gt_suffixerator'
env -i GT_MEM_BOOKKEEPING=on GTTESTDATA=${HOME}/gttestdata ./testsuite.rb \
       ${MC} -keywords 'gt_suffixerator and gttestdata' \
       -gttestdata ${GTTESTDATA}
env -i GT_MEM_BOOKKEEPING=on GTTESTDATA=${HOME}/gttestdata ./testsuite.rb \
       ${MC} -keywords gt_extractseq -gttestdata ${GTTESTDATA}
env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb ${MC} -keywords 'gt_trieins'
env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb ${MC} -keywords 'gt_ltrharvest'
env -i GT_MEM_BOOKKEEPING=on GTTESTDATA=${HOME}/gttestdata ./testsuite.rb \
       ${MC} -keywords 'gt_packedindex_at1MB' \
       -gttestdata ${GTTESTDATA}
env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb ${MC} -keywords 'gt_packedindex'
env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb ${MC} -keywords 'gt_greedyfwdmat'
env -i GT_MEM_BOOKKEEPING=on GTTESTDATA=${HOME}/gttestdata ./testsuite.rb \
       ${MC} -keywords 'gt_greedyfwdmat and gttestdata' \
       -gttestdata ${GTTESTDATA}
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
../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/Random-Small.fna
../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/Random.fna
../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/Atinsert.fna ../testdata/Random.fna
../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/TTT-small.fna
if test ! "X${GTTESTDATA}" = "X"
then
  AT=${GTTESTDATA}/Iowa/at1MB
  U8=${GTTESTDATA}/Iowa/U89959.fna
  ATK=${GTTESTDATA}/Iowa/at100K1
  GRUMBACH=${GTTESTDATA}/DNA-mix/Grumbach.fna
  ../scripts/rununique.sh ../bin/gt 10 20 ${U8} ${AT}
  for options in `${ALLOUTPUTOPTS}`
  do
    ../scripts/cmpdbfile.sh ${options} -pl -db ${ATK}
  done
  ../scripts/cmpdbfile.sh ${outoptions} -pl -db ${ATK} ${AT} ${GRUMBACH}/*.fna
fi
cd ..
