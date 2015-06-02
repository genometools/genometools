#!/bin/sh

set -e -x

len=0

if test $# -eq 1
then
  len=$1
else
  if test $# -ne 0
  then
    echo "Usage: $0 [len]"
    exit 1
  fi
fi


encseq2spm()
{
  gt encseq decode ${1} | head -n ${2} | grep -v '^>read_5' > tmp.fas
  gt encseq encode -indexname sfx tmp.fas 
  for threads in 1 4
  do
    GTCALL="env -i bin/gt -j ${threads} encseq2spm"
    ${GTCALL} -l 32 -singlescan 0 -ii sfx
    for len in 32 35
    do
      opts="-l ${len} -ii sfx -spm count -checksuftab"
      ${GTCALL} ${opts} 
    done
  done
}

found=0
for indexname in ${HOME}/seqcmpprojects/c22/c22-64bit /local/c22 ${HOME}/c22-64bit
do
  if test -f ${indexname}.esq
  then
    found=1
    break
  fi
done

if test $found -eq 0
then
  echo "$0: cannot find index for c22"
  exit 1
else
  if test $# -eq 0
  then
    for len in 100 300 600 900 15000 90000 300000
    do
      encseq2spm ${indexname} ${len}
    done
  else
    encseq2spm ${indexname} ${len}
  fi
fi
