#!/bin/sh

set -e -x

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

gt suffixerator -db $1  -dna -suf -tis -lcp -indexname sideeff
gt dev sfxmap -lcp -tis -suf -esa sideeff
