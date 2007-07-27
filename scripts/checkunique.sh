#!/bin/sh

set -e -x

if test $# -lt 4
then
  echo "Usage: $0 <min> <max> <queryfile> <file1> <file2> ..."
  exit 1
fi

rununique.sh ../bin/gt $* | grep -v '^#' > shit
rununique-test.sh $* | grep -v '^#' > shit-test
diff shit shit-test
