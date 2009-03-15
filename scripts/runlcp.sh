#!/bin/sh
checkerror() {
if test $? -ne 0
then
  echo "failure: ${cmd}"
  exit 1
else
  echo "okay: ${cmd}"
fi
}

set -e -x

for filename in `find testdata/ -name '*.fna'`
do
  valgrind.sh gt suffixerator -dna  -suf -lcp -tis -maxdepth -indexname sfx-idx -db ${filename}
done
