#!/bin/sh

set -e -x

for filename in `${GTDIR}/scripts/findfasta.rb`
do
  for t in 2 3 4 5 6
  do
    bin/gt -j $t suffixerator -parts 2 -dc 64 -suf -indexname sfx -db ${filename} 
    bin/gt dev sfxmap -esa sfx -suf
  done
done
