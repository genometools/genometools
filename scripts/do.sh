#!/bin/sh

for cmpcharbychar in 0 1
do
  for dir in fwd rev cpl rcl
  do
    extra=""
    if cmpcharbychar -eq 1
    then
      extra="-cmpcharbychar"
    fi
    cmd="gt suffixerator -dir ${dir} -parts 5 -dna -v -pl -bck -tis -suf -lcp -bwt -des -db ${AT} -showtime ${extra}"
    ${cmd}
    if test $? -ne 0
    then
      echo "failure: ${cmd}"
      exit 1
    else
      echo "success: ${cmd}"
    fi
  done
done
