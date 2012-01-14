#!/bin/sh

set -e -x

if test $# -ne 3
then
  echo "Usage: $0 <m|f> <minlen> <inputfile>"
  exit 1
fi

mirrored=$1
minlen=$2
inputfile=$3

env -i bin/gt sequniq -rev $inputfile > tmp.fas
sfxopt="-pl 2 -dna -lcp -suf -tis -ssp -db tmp.fas"
if test ${mirrored} == "m"
then
  idx="sfx-m"
else
  idx="sfx-f"
  env -i bin/gt encseq encode -indexname seq tmp.fas
  env -i bin/gt encseq decode -mirrored seq > tmp-m.fas
  mv tmp-m.fas tmp.fas
fi
sfxmapopt="-tis -suf -lcp -ssp -wholeleafcheck"
env -i bin/gt suffixerator -indexname ${idx}-$minlen -spmopt $minlen ${sfxopt}
env -i bin/gt dev sfxmap ${sfxmapopt} -esa ${idx}-$minlen
env -i bin/gt suffixerator -indexname ${idx} ${sfxopt}
env -i bin/gt dev sfxmap ${sfxmapopt} -esa ${idx}
env -i bin/gt repfind -spm -l $minlen -ii ${idx} > shit2
# env -i RUBYLIB=gtruby LD_LIBRARY_PATH=lib scripts/spmfilter.rb ${idx} shit1 > shit2
#env -i RUBYLIB=gtruby LD_LIBRARY_PATH=lib scripts/lcpintervals.rb -m $minlen ${idx} > shit3
#env -i RUBYLIB=gtruby LD_LIBRARY_PATH=lib scripts/lcpintervals.rb -m $minlen ${idx}-$minlen > shit4
env -i bin/gt dev sfxmap -stream-esq ${idx} storefirstcodes 32 > shit1
grep -v '^#' shit1 > shit5
#cmp -s shit2 shit3
#cmp -s shit3 shit4
cmp -s shit2 shit5 
rm -f shit[1-5]
