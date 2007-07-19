#!/bin/bash

if test $# -ne 1
then
  echo "Usage: $0 numofsequences"
  exit 1
fi

numofsequences=$1

INDEXNAME=../testdata/at1MB
TMPFILE=`mktemp TMP.XXXXXX` || exit 1
mkvtree.sh -indexname ${INDEXNAME} -db ${AT} -dna -pl -ois -tis
runVmatchprog.sh vseqselect -randomnum ${numofsequences} ${INDEXNAME} > ${TMPFILE}
splitmultifasta.pl TMP 60 0 ${TMPFILE}

filelist=`ls TMP-*`

runmerge.sh ${filelist}
rm -f ${TMPFILE} TMP-*
