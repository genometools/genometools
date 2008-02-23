#!/bin/sh

cmd="gt suffixerator -db ${AT} -dna -v -pl -bck -tis -suf -lcp -bwt -des"
${cmd}
if test $? -ne 0
then
  echo "failure: ${cmd}"
  exit 1
else
  echo "success: ${cmd}"
fi
