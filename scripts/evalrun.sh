#!/bin/sh
filenames=testdata/at1MB

suffixerator()
{
  printf "# PROGNAME ${filename} $*\n"
  ${RUNNER} gt suffixerator -showtime -indexname sfx-id -dna -tis -suf -db ${filename} $* | egrep 'PROGNAME|TIME overall'
} 

for filename in ${filenames}
do
  suffixerator -cmpcharbychar ""
  suffixerator ""
  suffixerator -cmpcharbychar -parts 3
  suffixerator -parts 3
  for dc in 8 16 32 64 128
  do
    suffixerator -cmpcharbychar -dc ${dc}
    suffixerator -dc ${dc}
  done
  rm -f sfx-idx.*
done
