#!/bin/bash

set -e -x

if test $# -lt 2
then
  echo "Usage: $0 <file1> <file2> ..."
  exit 1
fi

${GTDIR}/bin/gt suffixerator -indexname all -db $* -suf -lcp -pl 1
num=0
indexlist=""
for indexname in $*
do
  ${GTDIR}/bin/gt suffixerator -indexname midx${num} -db ${indexname} -suf -lcp -tis -pl 1
  indexlist="${indexlist} midx${num}"
  num=`expr ${num} + 1`
done
${GTDIR}/bin/gt dev mergeesa -indexname midx-all -ii ${indexlist}
cmp midx-all.suf all.suf
cmp midx-all.lcp all.lcp
cmp midx-all.llv all.llv
