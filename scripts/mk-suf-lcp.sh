#!/bin/sh

if test $# -eq 0
then
  echo "Usage: $0 <inputfile1> [inputfile2 ...]"
  exit 1
fi

env -i bin/gt suffixerator -tis -sds no -des no -md5 -suf -lcp -db $*
