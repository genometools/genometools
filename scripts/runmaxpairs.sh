#!/bin/sh

if test $# -lt 2
then
  echo "Usage: $0 minlength [filenames]"
  exit 1
fi

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

minlength=$1
shift
filenames=$*
VMMINI=vmatch-mini.x

for filename in ${filenames}
do
  echo "$0 ${minlength} ${filename}"
  checkerror "mkvtree.sh -db ${filename} -indexname mkvidx -dna -suf -tis -lcp -bwt -pl"
  checkerror "${VMMINI} ${minlength} mkvidx" > result.vm
  cleanhashlines result.vm
  checkerror "../bin/gt suffixerator -db ${filename} -indexname sfxidx -dna -suf -tis -lcp -pl"
  checkerror "../bin/gt dev maxpairs -l ${minlength} -ii sfxidx" > result.mp
  cleanhashlines result.mp
  cmp -s result.vm result.mp
done
