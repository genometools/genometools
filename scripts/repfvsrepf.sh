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
  sed -e '/^#/d' -e 's/[ ][ ]*/ /g' $1 > ${TMPFILE}
  mv ${TMPFILE} $1
}

sortlines()
{
  TMPFILE=`mktemp TMP.XXXXXX` || exit 1
  sort -n $1 > ${TMPFILE}
  mv ${TMPFILE} $1
}

extractlines()
{
  TMPFILE=`mktemp TMP.XXXXXX` || exit 1
  /sw/bin/gawk '/.*/ {print $1 " " $3 " " $4 " " $5 " " $7}' $1 > ${TMPFILE}
  mv ${TMPFILE} $1
}

if test $# -ne 2
then
  usage
  exit 1
fi

minlength=$1
filename=$2

GTDIR=/Users/stefan/genometools
checkerror "${GTDIR}/bin/gt suffixerator -db ${filename} -indexname sfxidx -dna -suf -tis -lcp -pl"
checkerror "${GTDIR}/bin/gt repfind -l ${minlength} -r -ii sfxidx" > result.gt
cleanhashlines result.gt
extractlines result.gt
sortlines result.gt
checkerror "/Users/stefan/bin-ops/i686-apple-darwin/repfind.x -allmax -l ${minlength} -r -noevalue -nodistance $filename" > result.rep
cleanhashlines result.rep
sortlines result.rep
checkerror "diff -w result.rep result.gt"
