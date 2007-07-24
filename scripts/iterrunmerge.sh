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
QUERY=${U8}
TMPFILE=`mktemp TMP.XXXXXX` || exit 1
mkdir -p ${INDEXDIR}
mkvtree.sh -indexname ${INDEXNAME} -db ${AT} -dna -pl -ois -tis
runVmatchprog.sh vseqselect -randomnum ${numofsequences} ${INDEXNAME} > ${TMPFILE}
splitmultifasta.pl TMP 60 0 ${TMPFILE}

filelist=`ls TMP-*`

./runmerge.sh ${U8} ${filelist}
rm -f ${TMPFILE} TMP-*
