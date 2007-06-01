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

checkerror "suffixerator.x -tis -indexname idx-sfx -pl 1 $*"
checkerror "mkvtree.sh -tis -indexname idx-mkv -dna -pl 1 -db $*"
checkerror "cmp -s idx-mkv.prj idx-sfx.prj"
