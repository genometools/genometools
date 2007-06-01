#!/bin/bash

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

for filename in $*
do
  for parts in 1 2 3
  do
    for sat in direct bit uchar ushort uint
    do
      checkerror "Checkmapped.sh -parts ${parts} -sat ${sat} -pl 3 ${filename}"
    done
    checkerror "suffixerator.x -tis -indexname sfx -parts 1 -sat uint64 -pl 3 ${filename}"
  done
done
