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

options="$*"

checkerror "${GTDIR}/bin/gt suffixerator -indexname /tmp/idx-sfx ${options}"
checkerror "mkvtree.sh -indexname /tmp/idx-mkv -dna ${options}"
checkerror "cmp -s /tmp/idx-mkv.prj /tmp/idx-sfx.prj"
