#!/usr/bin/env bash

set -e -x

if test $# -lt 4
then
  echo "Usage: $0 <min> <max> <queryfile> <file1> <file2> ..."
  exit 1
fi

INDEXDIR=../indexdir

mkdir -p ${INDEXDIR}

minvalue=$1
maxvalue=$2
queryfile=$3
shift
shift
shift
referencefiles=$*

numofreferencefiles=0
for referencefile in ${referencefiles}
do
  if test -f ${referencefile}
  then
    numofreferencefiles=`expr ${numofreferencefiles} + 1`
  else
    echo "$0: file \"${referencefile}\" does not exist"
    exit 1
  fi
done

if test ${numofreferencefiles} -eq 1
then
  ../bin/gt suffixerator -indexname ${INDEXDIR}/mkv-single -db ${referencefiles} -bwt  -pl 8
  ../bin/gt mkfmindex -noindexpos -fmout ${INDEXDIR}/fm-single -ii ${INDEXDIR}/mkv-single
  ../bin/gt suffixerator -indexname ${INDEXDIR}/fm-single -plain -smap ${INDEXDIR}/mkv-single.al1 -tis -pl 1 -db ${INDEXDIR}/mkv-single.bwt
  fmindexname=${INDEXDIR}/fm-single
else
  for referencefile in ${referencefiles}
  do
    ../bin/gt suffixerator -indexname ${INDEXDIR}/midx${num} -db ${referencefile} -suf -lcp -tis -pl 1
    indexlist="${indexlist} ${INDEXDIR}/midx${num}"
    num=`expr ${num} + 1`
  done
  ../bin/gt mkfmindex -noindexpos -fmout ${INDEXDIR}/fm-all -ii ${indexlist}
  ../bin/gt suffixerator -indexname ${INDEXDIR}/fm-all -plain -smap ${INDEXDIR}/fm-all.al1 -tis -pl 1 -db ${INDEXDIR}/fm-all.bwt
  fmindexname=${INDEXDIR}/fm-all
fi
../bin/gt uniquesub -fmi ${fmindexname} -query ${queryfile} -output sequence querypos -min ${minvalue} -max ${maxvalue}
