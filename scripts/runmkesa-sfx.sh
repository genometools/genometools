#!/bin/sh

set -e -x

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

time gt suffixerator -indexname sfx-idx -dna -v -suf -tis -showtime -pl -maxdepth -db $1 
time mkesa -p mkesa-idx -b D -g suf -v -d $1 
cmp -s sfx-idx.suf mkesa-idx.suf
