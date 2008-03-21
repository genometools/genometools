#!/bin/sh

for dir in fwd rev cpl rcl
do
  cmd="gt suffixerator -dir ${dir} -parts 5 -dna -v -pl -bck -tis -suf -lcp -bwt -des -db ${AT} -showtime"
${cmd}
  if test $? -ne 0
  then
    echo "failure: ${cmd}"
    exit 1
  else
    echo "success: ${cmd}"
  fi
done
