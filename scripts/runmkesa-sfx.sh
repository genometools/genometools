#!/bin/sh

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

# time mkesa -p mkesa-idx -b D -g suf -g lcp -v -d $1 
time gt suffixerator -dna -v -suf -tis -lcp -maxdepth -db $1 -showtime
