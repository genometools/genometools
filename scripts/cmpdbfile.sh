#!/bin/bash

if test $# -lt 1
then
  echo "Usage: $0 <files>"
  exit 1
fi

function checkerror() 
{
  $1
  if test $? -ne 0
  then
    echo "failure: ${1}"
    exit 1
  else
    echo "okay ${1}"
  fi
}

checkerror "${GTDIR}/bin/gt suffixerator -tis -indexname /tmp/idx-sfx -pl 1 -db $*"
checkerror "mkvtree.sh -tis -indexname /tmp/idx-mkv -dna -pl 1 -db $*"
TMPFILE1=`mktemp TMP.XXXXXX` || exit 1
grep -v 'integersize=' /tmp/idx-mkv.prj > ${TMPFILE1}
TMPFILE2=`mktemp TMP.XXXXXX` || exit 1
grep -v 'integersize=' /tmp/idx-sfx.prj > ${TMPFILE2}
checkerror "cmp -s ${TMPFILE1} ${TMPFILE2}"
