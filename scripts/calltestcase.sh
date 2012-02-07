#!/bin/sh 

if test $# -ne 1
then
  echo "Usage: $0 <keyword>"
  exit 1
fi

cd testsuite

if test $1 == "gt_encseq2spm"
then
  BK=off
else
  BK=on
fi
echo "run with GT_MEM_BOOKKEEPING=${BK}"
env -i GT_MEM_BOOKKEEPING=${BK}  ./testsuite.rb -gttestdata ${GTTESTDATA} \
         -keywords ${1}

cd ..
