#!/bin/sh
if test $# -ne 2
then
  echo "Usage: <file1> <file2>"
  exit 1
fi

set -e -x

WOTD="${TEACHDIR}/repertoire/exercises/effiziente_algorithmen_auf_sequenzen/suffix_trees/suffixtree_code/wotd.rb"
file1=$1
file2=$2

echo "# pairwise comparisons"
env -i ${GTDIR}/bin/gt suffixerator -db ${file1} -indexname idx1 -dna -suf -tis
env -i ${GTDIR}/bin/gt shulengthdist -ii idx1 -q ${file2} > shit.pairwise
env -i ${GTDIR}/bin/gt suffixerator -db ${file2} -indexname idx2 -dna -suf -tis
env -i ${GTDIR}/bin/gt shulengthdist -ii idx2 -q ${file1} >> shit.pairwise
env -i ${GTDIR}/bin/gt suffixerator -db ${file1} ${file2} -dna -suf -tis -lcp -indexname both 
echo "# multi comparisons"
env -i ${GTDIR}/bin/gt shulengthdist -ii both > shit.multi
diff shit.pairwise shit.multi
if test $? -ne 0
then
  echo "failure: cmp -s shit.pairwise shit.multi"
  exit 1
fi
# echo "`tail -1 tmp1`x`tail -1 tmp2`y" | ${WOTD} | dot -Tpdf > dot.pdf
# mkvtree -db tmp1 tmp2 -dna -suf -lcp -tis -indexname mkv-tmp
# vstree2tex -suf -lcp -tis -s mkv-tmp > tmp.tex
# latex tmp.tex
