#!/bin/sh

set -e -x

if test $# -ne 2
then
  echo "Usage: $0 minlength filename"
  exit 1
fi

minlength=$1
filename=$2


#mkvtree.sh -db ${filename} -indexname mkvidx -dna -suf -tis -lcp -bwt -pl
#/projects/vstree/src/vstree/src/Vmatch/vmatch.dbg.x -l ${minlength} mkvidx
bin/gt suffixerator -db ${filename} -indexname sfxidx -dna -suf -tis -lcp -pl
valgrind.sh bin/gt dev maxpairs -l ${minlength} -ii sfxidx
