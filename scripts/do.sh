#!/bin/sh
checkerror() {
if test $? -ne 0
then
  echo "failure: ${cmd}"
  exit 1
else
  echo "okay: ${cmd}"
fi
}


for cmpcharbychar in 0 1
do
  for dir in fwd rev cpl rcl
  do
    extra=""
    if cmpcharbychar -eq 1
    then
      extra="-cmpcharbychar"
    fi
    cmd="gt suffixerator -dir ${dir} -algbds 3 43 120 -parts 5 -dna -v -pl -bck -tis -suf -lcp -bwt -des -db ${AT} -showtime ${extra}"
    ${cmd}
    checkerror
    cmd="gt dev sfxmap -tis -suf -lcp -bwt -des -esa at1MB"
    ${cmd}
    checkerror
  done
done
