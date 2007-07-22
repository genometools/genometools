#!/bin/sh

set -e -x

if test $# -lt 3
then
  echo "Usage: $0 <queryfile> <file1> <file2> ..."
  exit 1
fi

PREFIX=../testdata

queryfile=$1
shift
referencefiles=$*
../bin/gt suffixerator -indexname ${PREFIX}/all -db ${referencefiles} -suf -lcp -pl 1
num=0
indexlist=""
for filename in ${referencefiles}
do
  ../bin/gt suffixerator -indexname ${PREFIX}/midx${num} -db ${filename} -suf -lcp -tis -pl 1
  indexlist="${indexlist} ${PREFIX}/midx${num}"
  num=`expr ${num} + 1`
done
../bin/gt dev mergeesa -indexname ${PREFIX}/midx-all -ii ${indexlist}
cmp ${PREFIX}/midx-all.suf ${PREFIX}/all.suf
cmp ${PREFIX}/midx-all.lcp ${PREFIX}/all.lcp
cmp ${PREFIX}/midx-all.llv ${PREFIX}/all.llv
../bin/gt mkfmindex -noindexpos -fmout ${PREFIX}/fm-all -ii ${indexlist}
../bin/gt suffixerator -indexname ${PREFIX}/fm-all -plain -smap ${PREFIX}/fm-all.al1 -tis -pl 1 -db ${PREFIX}/fm-all.bwt
../bin/gt uniquesub -fmi ${PREFIX}/fm-all -query ${queryfile} -output sequence querypos -min 10 -max 10
