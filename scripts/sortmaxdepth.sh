#!/bin/sh

runall()
{
  bin/gt suffixerator -tis -db ${1} -indexname esa
  for maxinsertionsort in 1 2 5 10
  do
    for maxbltriesort in 10 20 30 100
    do
      for maxdepth in 2 3 4 5 6 7 8 9 10 11 12 13
      do
        bin/gt dev sfxmap -sortmaxdepth ${maxdepth} -esa esa -v \
                          -algbds ${maxinsertionsort} ${maxbltriesort} \
                                  1000
      done
      bin/gt suffixerator -db ${1} -indexname esa2 -suf -v \
                          -algbds ${maxinsertionsort} ${maxbltriesort} \
                                  1000
      bin/gt dev sfxmap -suf -esa esa2 
    done
  done
}

set -e -x

for filename in `ls testdata/[A-Z]*.fna`
do
  if test ${filename} != "testdata/corruptpatternfile.fna"
  then
    runall ${filename}
  fi
done
