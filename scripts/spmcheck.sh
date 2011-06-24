#!/bin/sh

set -e -x

if test $# -ne 2
then
  echo "Usage: $0 <minlen> <inputfile>"
  exit 1
fi

sfxopt="-dna -lcp -suf -tis -ssp -db $2"
sfxmapopt="-tis -suf -lcp -ssp -wholeleafcheck"
env -i bin/gt suffixerator -mirrored -indexname sfx-$1 -spmopt $1 ${sfxopt}
env -i bin/gt dev sfxmap ${sfxmapopt} -esa sfx-$1
env -i bin/gt suffixerator -mirrored -indexname sfx ${sfxopt}
env -i bin/gt dev sfxmap ${sfxmapopt} -esa sfx
env -i bin/gt repfind -l $1 -ii sfx > shit1
scripts/spmfilter.rb sfx shit1 > shit2
scripts/lcpintervals.rb -m $1 sfx > shit3
scripts/lcpintervals.rb -m $1 sfx-$1 > shit4
cmp -s shit2 shit3
cmp -s shit3 shit4
rm -f shit[1-4]
