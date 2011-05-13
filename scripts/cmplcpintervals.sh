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
  cmd="scripts/lcpintervals.rb sfx"
  ${cmd} > tmp1.result
  checkerror
  cmd="env -i bin/gt dev sfxmap -enumlcpintervals -esa sfx"
  ${cmd} > tmp2.result
  checkerror
  cmd="cmp -s tmp1.result tmp2.result"
  ${cmd}
  checkerror
done
rm -f tmp1.result tmp2.result
