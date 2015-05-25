#!/bin/sh

set -e -x

for filename in `${GTDIR}/scripts/findfasta.rb`
do
  gt suffixerator -db ${filename} -suftabuint -lcp -suf -indexname sfx -compressedoutput
  gt dev sfxmap -esa sfx -compresslcp
  gt dev sfxmap -esa sfx -compressedesa
done
