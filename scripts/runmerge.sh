#!/bin/sh

set -e -x

if test $# -lt 3
then
  echo "Usage: $0 <queryfile> <file1> <file2> ..."
  exit 1
fi

INDEXDIR=../indexdir
mkdir -p ${INDEXDIR}

queryfile=$1
shift
referencefiles=$*
../bin/gt suffixerator -dna -algbds 3 40 120 -indexname ${INDEXDIR}/all -db ${referencefiles} -suf -lcp -tis -pl
num=0
indexlist=""
for filename in ${referencefiles}
do
  ../bin/gt suffixerator -dna -v -cmpcharbychar -indexname ${INDEXDIR}/midx${num} -db ${filename} -suf -lcp -tis -pl
  ../bin/gt dev sfxmap -tis -suf -lcp  ${INDEXDIR}/midx${num}
  indexlist="${indexlist} ${INDEXDIR}/midx${num}"
  num=`expr ${num} + 1`
done
../bin/gt dev mergeesa -indexname ${INDEXDIR}/midx-all -ii ${indexlist}
cmp ${INDEXDIR}/midx-all.suf ${INDEXDIR}/all.suf
cmp ${INDEXDIR}/midx-all.lcp ${INDEXDIR}/all.lcp
cmp ${INDEXDIR}/midx-all.llv ${INDEXDIR}/all.llv
../bin/gt mkfmindex -noindexpos -fmout ${INDEXDIR}/fm-all -ii ${indexlist}
../bin/gt suffixerator -algbds 3 40 120 -indexname ${INDEXDIR}/fm-all -plain -smap ${INDEXDIR}/fm-all.al1 -tis -pl 1 -db ${INDEXDIR}/fm-all.bwt
../bin/gt uniquesub -fmi ${INDEXDIR}/fm-all -query ${queryfile} -output sequence querypos -min 10 -max 10
