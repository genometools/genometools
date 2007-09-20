#!/bin/sh

if test $# -lt 3
then
  echo "Usage: $0 <suffixerator-options>"
  exit 1
fi

checkerror() 
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

#VALGRIND=valgrind.sh

checksfxmap()
{
  cmd="${VALGRIND} ../bin/gt dev sfxmap -v $1"
  TMPFILE=`mktemp TMP.XXXXXX` || exit 1
  ${cmd} > ${TMPFILE}
  if test $? -ne 0
  then
    echo "failure: ${cmd}"
    exit 1
  else
    echo "success: ${cmd}"
  fi
  grep -v '^#' ${TMPFILE} > $1-sfx.prj
  rm -f ${TMPFILE}
  checkerror "cmp -s $1.prj $1-sfx.prj"
}

checkerror "${VALGRIND} ../bin/gt suffixerator -tis -des -suf -bwt -lcp -indexname /tmp/sfxidx $*"
checkerror "checksfxmap /tmp/sfxidx"
checkerror "../bin/gt dev sfxmap -v -stream /tmp/sfxidx"
