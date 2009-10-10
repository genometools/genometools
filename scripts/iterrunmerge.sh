#!/bin/sh

set -e -x

if test $# -ne 1
then
  echo "Usage: $0 numofsequences"
  exit 1
fi

numofsequences=$1

if test ! "X${GTTESTDATA}" = "X"
then
  QUERY=../testdata/U89959_genomic.fas
  AT=../testdata/at1MB
  TMPFILE=`mktemp TMP.XXXXXX` || exit 1
  ../scripts/seqselect.rb ${numofsequences} ${AT} > ${TMPFILE}
  ../scripts/splitmultifasta.rb TMP 0 ${TMPFILE}
  filelist=`ls TMP-*`
  ../scripts/runmerge.sh ${QUERY} ${filelist}
  rm -f ${TMPFILE} TMP-*
fi
