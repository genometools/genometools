#!/bin/sh
if test $# -eq 0
then
  echo "Usage: $0 <command to debug>"
  exit 1
fi
TMPFILE=`mktemp -t tmpfile.XXXXXX`
shift
echo "break exit" >> ${TMPFILE}
echo "run $*" >> ${TMPFILE}
echo "bt" >> ${TMPFILE}
cat ${TMPFILE} | gdb ${GTDIR}/bin/gt
rm -f ${TMPFILE}
