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
  cmd="env -i bin/gt suffixerator -tis -suf -lcp -dna -indexname sfx -db ${filename}"
  ${cmd}
  checkerror
  cmd="env -i RUBYLIB=gtruby LD_LIBRARY_PATH=lib scripts/lcpintervals.rb --itv sfx"
  ${cmd} > itvs.result1
  checkerror
  cmd="env -i bin/gt dev sfxmap -enumlcpitvs -esa sfx"
  ${cmd} > itvs.result2
  checkerror
  cmd="cmp -s itvs.result1 itvs.result2"
  ${cmd}
  checkerror
  cmd="env -i RUBYLIB=gtruby LD_LIBRARY_PATH=lib scripts/lcpintervals.rb --tree sfx"
  ${cmd} > itvtree.result1
  checkerror
  cmd="env -i bin/gt dev sfxmap -enumlcpitvtree -esa sfx"
  ${cmd} > itvtree.result2
  checkerror
  cmd="env -i bin/gt dev sfxmap -enumlcpitvtreeBU -esa sfx"
  ${cmd} > itvtree.result3
  checkerror
  cmd="cmp -s itvtree.result1 itvtree.result2"
  ${cmd}
  checkerror
  cmd="cmp -s itvtree.result2 itvtree.result3"
  ${cmd}
  checkerror
done
rm -f itvs.result1 itvs.result2
rm -f itvtree.result1 itvtree.result2 itvtree.result3
