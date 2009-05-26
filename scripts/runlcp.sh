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

suffixeratornoidx()
{
  ${RUNNER} gt suffixerator -v -showtime -dna -tis -lcp -suf -des -ssp -db ${filename} $*
}

suffixerator()
{
  ${RUNNER} gt suffixerator -v -showtime -dna -tis -lcp -suf -des -ssp -db ${filename} -indexname sfx-idx $*
}

suffixeratoronlysuf()
{
  ${RUNNER} gt suffixerator -v -showtime -dna -tis -suf -db ${filename} -indexname sfx-idx $*
}

sfxmap()
{
  gt dev sfxmap -lcp -suf $*
}

sfxmaponlysuf()
{
  gt dev sfxmap -suf $*
}

for filename in ${filenames}
do
  suffixerator ""
  sfxmap sfx-idx
  suffixerator -dir rev
  sfxmap sfx-idx
  suffixerator -maxdepth
  sfxmap sfx-idx
  maxdepth=`grep '^prefixlength=' sfx-idx.prj | sed -e 's/prefixlength=//'`
  maxdepth=`expr ${maxdepth} \* 2`
  suffixeratornoidx -maxdepth ${maxdepth} -indexname sfx-idx${maxdepth}
  sfxmap sfx-idx${maxdepth}
  suffixerator -parts 3
  sfxmap sfx-idx
  suffixerator -parts 3 -maxdepth
  sfxmap sfx-idx
  suffixerator -parts 3 -maxdepth he
  sfxmap sfx-idx
  suffixerator -parts 3 -maxdepth abs
  sfxmap sfx-idx
  suffixeratoronlysuf -algbds 20 20 100 -cmpcharbychar -dc 8
  sfxmaponlysuf sfx-idx
  suffixeratoronlysuf -algbds 20 20 100 -cmpcharbychar -dc 32
  sfxmaponlysuf sfx-idx
  suffixeratoronlysuf -algbds 20 20 100 -dc 32
  sfxmaponlysuf sfx-idx
  suffixeratoronlysuf -dc 32
  sfxmaponlysuf sfx-idx
  suffixeratoronlysuf -algbds 20 20 100 -dc 8
  sfxmaponlysuf sfx-idx
  suffixeratoronlysuf -cmpcharbychar -dc 16
  sfxmaponlysuf sfx-idx
  suffixeratoronlysuf -cmpcharbychar -dc 8
  sfxmaponlysuf sfx-idx
  ${RUNNER} gt suffixerator -v -showtime -smap Transab -tis -suf -dc 64 -db testdata/fib25.fas.gz
  rm -f sfx-idx.* sfx-idx${maxdepth}.*
done
echo "${filenames}"
