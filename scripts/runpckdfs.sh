#!/bin/sh

set -e -x

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi
inputfile=$1

gt suffixerator -db ${inputfile} -dna -suf -indexname esa
gt suffixerator -db ${inputfile} -dna -suf -lcp -dir rev -indexname esa-rev
gt packedindex mkindex -db ${inputfile} -dna -indexname pck
gt dev sfxmap -esa esa -pck pck -cmpsuf -v
gt dev sfxmap -esa esa-rev -pck pck -cmplcp -v
