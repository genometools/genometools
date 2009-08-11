#!/bin/sh

if test $# -ne 2
then
  echo "Usage: $0 <fromcasenum> <tocasenum>"
  exit 1
fi

apply()
{
  echo $1
}


testnum=$1
while test ${testnum} -le $2
do
  filepath="testsuite/stest_testsuite/test${testnum}"
  apply ${filepath}
  testnum=`expr ${testnum} + 1`
done
