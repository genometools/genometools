#!/bin/sh 

checkerror() {
if test $? -ne 0
then
  echo "failure: ${cmd}"
  exit 1
fi
}

if test $# -eq 1
then
  filenamelist="$1"
else
  filenamelist=`find testdata/ -name '[A-Z]*.fna'`
fi

for filename in ${filenamelist}
do
  echo "${filename}"
  cmd="env -i bin/gt suffixerator -db ${filename} -tis -suf -lcp -dna -indexname sfx"
  ${cmd}
  checkerror
  cmd="scripts/lcpintervals.rb itv sfx"
  ${cmd} > itvs.result1
  checkerror
  cmd="env -i bin/gt dev sfxmap -enumlcpintervals -esa sfx"
  ${cmd} > itvs.result2
  checkerror
  cmd="cmp -s itvs.result1 itvs.result2"
  ${cmd}
  checkerror
  cmd="scripts/lcpintervals.rb debugtree sfx"
  ${cmd} > tmp.result
  checkerror
  grep -v '^#' tmp.result > itvtree.result1
  cmd="env -i bin/gt dev sfxmap -enumlcpintervaltree -esa sfx"
  ${cmd} > itvtree.result2
  checkerror
  cmd="cmp -s itvtree.result1 itvtree.result2"
  ${cmd}
  checkerror
done
rm -f itvs.result1 itvs.result2
rm -f itvtree.result1 itvtree.result2
