#!/bin/sh 

if test $# -ne 1
then
  echo "Usage: $0 <keyword>"
  exit 1
fi

cd testsuite

./testsuite.rb -gttestdata ${GTTESTDATA} -keywords ${1}

cd ..
