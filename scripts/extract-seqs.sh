#!/bin/sh

if test $# -lt 2
then
  echo "Usage: $0 <fastafile> seqnum1 [seqnum2 ...]"
  exit 1
fi

inputfile=$1
shift
for seqnum in $*
do
  env -i bin/gt seq -width 70 -showseqnum `expr ${seqnum} + 1` ${inputfile}
done
