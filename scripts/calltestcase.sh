#!/bin/sh 

if test $# -eq 0
then
  echo "Usage: $0 <keyword1> [keyword2 ...]"
  exit 1
fi

cd testsuite

for keyword in $*
do
  echo "run testsuite for keyword ${keyword}"
  env -i ./testsuite.rb -gttestdata ${GTTESTDATA} -keywords $keyword
done

cd ..
