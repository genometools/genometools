#!/bin/sh
set -e -x

if test $# -eq 0
then
  filenames=`find testdata/ -name '*.fna'`
else
  if test $1 == 'valgrind'
  then
    VALGRIND=valgrind.sh
    shift
  fi
  if test $# -eq 0
  then
    filenames=`find testdata/ -name '*.fna'`
  else
    filenames=$*
  fi
fi


for filename in ${filenames}
do
  ${VALGRIND} gt suffixerator -v -lcp -tis -suf -des -ssp -maxdepth -indexname sfx-idx -db ${filename}
  gt dev sfxmap -suf sfx-idx
  maxdepth=`grep '^prefixlength=' sfx-idx.prj | sed -e 's/prefixlength=//'`
  maxdepth=`expr ${maxdepth} \* 2`
  ${VALGRIND} gt suffixerator -v -lcp -tis -suf -des -ssp -maxdepth ${maxdepth} -indexname sfx-idx${maxdepth} -db ${filename}
  cmp sfx-idx${maxdepth}.suf sfx-idx.suf 
  gt dev sfxmap -suf sfx-idx${maxdepth}
done
