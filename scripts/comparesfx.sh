#!/bin/sh

for filename in `scripts/enumseqfiles.sh`
do
  # echo "${filename}"
  bin/gt suffixerator -db ${filename} -tis -v > tmp1
  /local/kurtz/tallymer-test/bin/gt suffixerator -algbds 3 43 120 -db ${filename} -tis -v > tmp2
  diff tmp1 tmp2
done
