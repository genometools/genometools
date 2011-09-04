#!/bin/sh 

if test $# -ne 1
then
  echo "Usage: $0 <keyword>"
  exit 1
fi

cd testsuite

env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb -gttestdata ${GTTESTDATA} \
       -keywords gt_${1}

cd ..
