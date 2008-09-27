#!/bin/sh

set -e -x

cerr()
{
  $*
  if [ $? -ne 0 ]
  then
    echo "failure: $*"
    exit 1
  fi
}

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

mersize=7
inputfile=$1

cerr "bin/gt suffixerator -db ${inputfile} -tis -suf -lcp -pl -dna -indexname sfxidx"
cerr "env -i bin/gt tallymer mkindex -mersize 7 sfxidx" > tmp2
cerr "mkvtree.x -db ${inputfile} -tis -suf -lcp -pl -dna -indexname mkvidx"
cerr "tallymer-mkindex -mersize 7 mkvidx" | grep -v '^#' > tmp1
cerr "cmp -s tmp1 tmp2"
