#!/bin/sh

set -e -x

if test $# -ne 1
then
  echo "Usage: $0 <filename>"
  exit 1
fi

runsfx()
{
  env -i bin/gt suffixerator -dna -suf -lcp -tis -spmopt 20 $*
}

cp $1 both.fna
reversecompl.pl 60 $1 >> both.fna

runsfx -indexname sfx1 -db $1 -mirrored
runsfx -indexname sfx2 -db both.fna
cmp -s sfx1.suf sfx2.suf
