#!/bin/sh

#set -e -x

usage()
{
  echo "Usage: $0 minlength [-q query] file1 [file2 ...]"
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
  sed -e '/^#/d' $1 > ${TMPFILE}
  mv ${TMPFILE} $1
}

if test $# -lt 2
then
  usage
  exit 1
fi

queryfile=""
minlength=$1
if test "$2" = "-q"
then
  if test $# -lt 4
  then
    usage
    exit 1
  fi
  queryfile=$3 
  shift
  shift
fi
shift
filenames=$*

for filename in ${filenames}
do
  echo "$0 ${minlength} ${filename}"
  checkerror "../bin/gt suffixerator -algbds 3 43 120 -db ${filename} -indexname sfxidx -dna -suf -tis -lcp -pl"
  if test "X${queryfile}" = "X"
  then
    checkerror "../bin/gt repfind -l ${minlength} -ii sfxidx" > result.mp
    cleanhashlines result.mp
  else
    checkerror "../bin/gt repfind -l ${minlength} -q ${queryfile} -ii sfxidx" > result.mp
    cleanhashlines result.mp
  fi
  resultfile="${GTTESTDATA}/repfind-result/`basename ${filename}`.result"
  checkerror "cmp -s ${resultfile} result.mp"
done
