#!/bin/bash

if test $# -ne 1
then
  echo "Usage: $0 numofsequences"
  exit 1
fi

numofsequences=$1

TMPFILE=`mktemp TMP.XXXXXX` || exit 1
mkvtree.sh -db ${AT} -dna -pl -ois -tis
runVmatchprog.sh vseqselect -randomnum ${numofsequences} at1MB > ${TMPFILE}
splitmultifasta.pl TMP 60 0 ${TMPFILE}

filelist=`ls TMP-*`

runmerge.sh ${filelist}
rm -f ${TMPFILE} TMP-*
