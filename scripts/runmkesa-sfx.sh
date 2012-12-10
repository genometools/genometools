#!/bin/sh

# set -e -x
TMP=.time.tmp

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

/usr/bin/time -o ${TMP} -f "# TIME overall %U" mkesa -p mkesa-idx -b D -g suf -v -d $1 |\
                             mkesa-fmt.rb
cat ${TMP}
rm -f ${TMP}
# mkesa -p mkesa-idx -b D -g suf -v -d $1
# time gt suffixerator -indexname sfx-idx -dna -v -suf -tis -showtime -pl -dc 64 -db $1 
#cmp -s sfx-idx.suf mkesa-idx.suf
