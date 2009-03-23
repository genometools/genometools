#!/bin/sh
set -e -x

if test $# -eq 0
then
  filenames="`find testdata/ -name '*.fna'` testdata/at1MB"
else
  if test $1 == 'valgrind'
  then
    RUNNER=valgrind.sh
    shift
  else
    RUNNER=time
  fi
  if test $# -eq 0
  then
    filenames="`find testdata/ -name '*.fna'` testdata/at1MB"
  else
    filenames=$*
  fi
fi

suffixerator()
{
  ${RUNNER} gt suffixerator -lcp -tis -suf -des -ssp -db ${filename} $*
}

sfxmap()
{
  gt dev sfxmap -lcp -suf $*
}

for filename in ${filenames}
do
  suffixerator -indexname sfx-idx 
  sfxmap sfx-idx
  suffixerator -maxdepth -indexname sfx-idx 
  sfxmap sfx-idx
  maxdepth=`grep '^prefixlength=' sfx-idx.prj | sed -e 's/prefixlength=//'`
  maxdepth=`expr ${maxdepth} \* 2`
  suffixerator -maxdepth ${maxdepth} -indexname sfx-idx${maxdepth}
  sfxmap sfx-idx${maxdepth}
  rm -f sfx-idx.* sfx-idx${maxdepth}.*
done
echo "${filenames}"
