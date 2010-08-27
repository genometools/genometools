#!/bin/sh
if test $# -ne 2
then
  echo "Usage: <file1> <file2>"
  exit 1
fi

WOTD="${TEACHDIR}/repertoire/exercises/effiziente_algorithmen_auf_sequenzen/suffix_trees/suffixtree_code/wotd.rb"
file1=$1
file2=$2

echo "# pairwise comparisons"
env -i bin/gt suffixerator -db ${file1} -dna -suf -tis
printf "0 1 " > shit.pairwise
env -i bin/gt shulengthdist -ii ${file1} -q ${file2} >> shit.pairwise
env -i bin/gt suffixerator -db ${file2} -dna -suf -tis
printf "1 0 " >> shit.pairwise
env -i bin/gt shulengthdist -ii ${file2} -q ${file1} >> shit.pairwise
env -i bin/gt suffixerator -db ${file1} ${file2} -dna -suf -tis -lcp -indexname both 
echo "# multi comparisons"
env -i bin/gt shulengthdist -ii both > shit.multi
cmp -s shit.pairwise shit.multi
# echo "`tail -1 tmp1`x`tail -1 tmp2`" | ${WOTD} | dot -Tpdf > dot.pdf
