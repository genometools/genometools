#!/bin/sh

if test $# -lt 1
then
  echo "Usage: $0 <files>"
  exit 1
fi

checkerror() 
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

comparefiles()
{
  TMPFILE1=`mktemp ../testsuite/TMP.XXXXXX` || exit 1
  egrep -v 'prefixlength|integersize|readmode|realspecialranges' $1 > ${TMPFILE1}
  TMPFILE2=`mktemp ../testsuite/TMP.XXXXXX` || exit 1
  egrep -v 'prefixlength|integersize|readmode|realspecialranges' $2 > ${TMPFILE2}
  checkerror "cmp -s ${TMPFILE1} ${TMPFILE2}"
  rm -f ${TMPFILE1} ${TMPFILE2}
}

options="$*"

checkerror "../bin/gt suffixerator -indexname /tmp/idx-sfx -des ${options}"
checkerror "mkvtree.sh -indexname /tmp/idx-mkv -dna ${options}"
comparefiles /tmp/idx-mkv.prj /tmp/idx-sfx.prj
