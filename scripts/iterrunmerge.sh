#!/bin/sh

set -e -x

if test $# -ne 1
then
  echo "Usage: $0 numofsequences"
  exit 1
fi

numofsequences=$1

INDEXDIR=../indexdir
INDEXNAME=${INDEXDIR}/at1MB
if test ! "X${GTTESTDATA}" = "X"
then
  QUERY=${GTTESTDATA}/Iowa/U89959.fna
  AT=${GTTESTDATA}/Iowa/at1MB
  TMPFILE=`mktemp TMP.XXXXXX` || exit 1
  mkdir -p ${INDEXDIR}
  mkvtree.sh -indexname ${INDEXNAME} -db ${AT} -dna -pl -ois -tis
  runVmatchprog.sh vseqselect -randomnum ${numofsequences} ${INDEXNAME} > ${TMPFILE}
  splitmultifasta.pl TMP 60 0 ${TMPFILE}
  filelist=`ls TMP-*`
  ../scripts/runmerge.sh ${QUERY} ${filelist}
  rm -f ${TMPFILE} TMP-*
fi
