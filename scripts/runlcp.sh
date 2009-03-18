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
  ${VALGRIND} gt suffixerator -tis -suf -des -ssp -lcp -bwt -bck  -maxdepth -indexname sfx-idx -db ${filename}
  gt dev sfxmap -suf -lcp sfx-idx
done
