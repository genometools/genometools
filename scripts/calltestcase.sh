#!/bin/sh 

cd testsuite

if test "${TMPDIR}" == ""
then
  TMPDIR="."
fi

TMPFILE=`mktemp ${TMPDIR}/TMP.XXXXXX` || exit 1
if test $# -eq 0
then
  echo ${GTTESTKEYWORDS} | tr ' ' '\n' > ${TMPFILE}
else
  touch ${TMPFILE}
  for keyword in $*
  do
    echo $keyword >> ${TMPFILE}
  done
fi

for keyword in `cat ${TMPFILE}`
do
  grep -q $keyword *.rb
  if test $? -ne 0
  then
    echo "$0: keyword ${keyword} not in testsuite/*.rb"
    exit 1
  fi
done

for keyword in `cat ${TMPFILE}`
do
  echo "run testsuite for keyword ${keyword}"
  env -i ./testsuite.rb -gttestdata ${GTTESTDATA} -keywords $keyword
done

rm -f ${TMPFILE}

cd ..
