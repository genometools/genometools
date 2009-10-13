#!/bin/sh

#set -e -x

usage()
{
  echo "Usage: $0 minlength inputfile"
}

checkerror() 
{
  $1
  if test $? -ne 0
  then
    echo "failure: ${1}"
    exit 1
  fi
}

cleanhashlines()
{
  TMPFILE=`mktemp TMP.XXXXXX` || exit 1
  sort $1 | sed -e '/^#/d' > ${TMPFILE}
  mv ${TMPFILE} $1
}

if test $# -ne 2
then
  usage
  exit 1
fi

minlength=$1
filename=$2

checkerror "${GTDIR}/bin/gt suffixerator -db ${filename} -indexname sfxidx -dna -suf -tis -lcp -pl"
checkerror "valgrind.sh ${GTDIR}/bin/gt repfind -l ${minlength} -r -ii sfxidx" > result.gt
cleanhashlines result.gt
checkerror "repfind.x -l ${minlength} -r -noevalue -nodistance $filename" > result.rep
cleanhashlines result.rep
checkerror "diff result.rep result.gt"
