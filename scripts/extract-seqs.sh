#!/bin/sh

if test $# -lt 2
then
  echo "Usage: $0 <fastafile|encseq> seqnum1 [seqnum2 ...]"
  exit 1
fi

inputfile=$1
shift
for seqnum in $*
do
  if test -f ${inputfile}.esq
  then
    env -i bin/gt encseq decode -seq ${seqnum} ${inputfile}
  else
    env -i bin/gt seq -width 70 -showseqnum `expr ${seqnum} + 1` ${inputfile}
  fi
done
