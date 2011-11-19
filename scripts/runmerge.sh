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
SFXOPTS="-dna -algbds 3 43 120 -suf -lcp -tis -pl"
../bin/gt suffixerator ${SFXOPTS} -indexname ${INDEXDIR}/all -db ${referencefiles}
num=0
indexlist=""
for filename in ${referencefiles}
do
  ../bin/gt suffixerator ${SFXOPTS} -indexname ${INDEXDIR}/midx${num} -db ${filename}
  ../bin/gt dev sfxmap -tis -suf -lcp -esa ${INDEXDIR}/midx${num}
  indexlist="${indexlist} ${INDEXDIR}/midx${num}"
  num=`expr ${num} + 1`
done
../bin/gt dev mergeesa -indexname ${INDEXDIR}/midx-all -ii ${indexlist}
cmp ${INDEXDIR}/midx-all.suf ${INDEXDIR}/all.suf
cmp ${INDEXDIR}/midx-all.lcp ${INDEXDIR}/all.lcp
cmp ${INDEXDIR}/midx-all.llv ${INDEXDIR}/all.llv
../bin/gt mkfmindex -noindexpos -fmout ${INDEXDIR}/fm-all -ii ${indexlist}
../bin/gt suffixerator -indexname ${INDEXDIR}/fm-all -plain -des no -ssp no -sds no -smap ${INDEXDIR}/fm-all.al1 -tis -db ${INDEXDIR}/fm-all.bwt
../bin/gt uniquesub -fmi ${INDEXDIR}/fm-all -query ${queryfile} -output sequence querypos -min 10 -max 10
