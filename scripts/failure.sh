#!/bin/sh

set -e -x

if test $# -ne 2
then
  echo "Usage: $0 <Referencefile> <Queryfile>"
  exit 1
fi  

reference=$1
query=$2

# exec in gt root dir

# mkvtree.x -db ${reference} -indexname esa-mkv -dna -pl -lcp -suf -tis 
# vstree2tex.x -tis -suf -s esa-mkv > tmp.tex
gt suffixerator -tis -suf -dna -pl -db ${reference} -indexname esa-fwd
gt tagerator -e 2 -maxocc 20 -nospecials -rw -esa esa-fwd -q ${query}
gt packedindex mkindex -tis -indexname pck-rev -db ${reference} -sprank -dna -pl -bsize 10 -locfreq 32 -dir rev
gt tagerator -e 2 -maxocc 20 -nospecials -rw -pck pck-rev -q ${query}
