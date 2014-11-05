#!/bin/sh -e

if [ $# -eq 3 ]
then
  bn=$(basename $3)
  $GTDIR/bin/gt seqfilter -minlength $1 $3 | \
    $GTDIR/bin/gt shredder -sample 0.1 \
    -clipdesc \
    -minlength $1 \
    -maxlength $1 \
    -coverage $2 | \
    $GTDIR/bin/gt seqfilter -minlength $1 -width 70 > \
    ${bn%%.fas}_queries_$1_${2}x.fas
  exit 0
fi

echo "Usage: " $0 " SHREDDLEN COVERAGE INFILE"
exit 1

