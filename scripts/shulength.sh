#!/bin/sh
if test $# -ne 2
then
  echo "Usage: <file1> <file2>"
  exit 1
fi

file1=$1
file2=$2

echo "pairwise comparisons"
env -i bin/gt suffixerator -db ${file1} -dna -suf -tis
printf "0 1 "
env -i bin/gt shulengthdist -ii ${file1} -q ${file2}
env -i bin/gt suffixerator -db ${file2} -dna -suf -tis
printf "1 0 "
env -i bin/gt shulengthdist -ii ${file2} -q ${file1}
env -i bin/gt suffixerator -db ${file1} ${file2} -dna -suf -tis -lcp -indexname both
echo "multi comparisons"
env -i bin/gt shulengthdist -ii both
