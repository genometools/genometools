#!/bin/bash

if test $# -lt 3
then
  echo "Usage: $0 <suffixerator-options>"
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

function checksfxmap()
{
  cmd="sfxmap.x $1"
  TMPFILE=`mktemp TMP.XXXXXX` || exit 1
  ${cmd} > ${TMPFILE}
  if test $? -ne 0
  then
    echo "failure: sfxmap.x ${1}"
    exit 1
  fi
  grep -v '^#' ${TMPFILE} > $1-sfx.prj
  rm -f ${TMPFILE}
  checkerror "cmp -s $1.prj $1-sfx.prj"
}

checkerror "suffixerator.x -tis -suf -indexname sfx $*"
checkerror "checksfxmap sfx"
