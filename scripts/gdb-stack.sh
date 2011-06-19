#!/bin/sh

TMPFILE=`mktemp -t tmpfile.XXXXXX`
shift
echo "break exit" >> ${TMPFILE}
echo "run $*" >> ${TMPFILE}
echo "bt" >> ${TMPFILE}
cat ${TMPFILE} | gdb bin/gt
rm -f ${TMPFILE}
