#!/bin/sh
set -e -x

if test $# -eq 0
then
  filenames=`find testdata/ -name '*.fna'`
else
  filenames=$*
fi

for filename in ${filenames}
do
  gt suffixerator -dna -v -tis -suf -des -ssp -lcp -bwt -bck  -maxdepth -indexname sfx-idx -db ${filename}
  gt dev sfxmap -suf -lcp sfx-idx
done
