#!/bin/sh

cmd="gt suffixerator -db ${AT} -parts 5 -dna -v -pl -bck -tis -suf -lcp -bwt -des"
${cmd}
if test $? -ne 0
then
  echo "failure: ${cmd}"
  exit 1
else
  echo "success: ${cmd}"
fi
