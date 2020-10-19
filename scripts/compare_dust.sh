#!/bin/sh

set -e -x

which -s dustmasker
if test $? -ne 0
then
  echo "$0: cannot find dustmasker"
  exit 1
fi
DUSTMASKER=dustmasker

dust_compare()
{
  inputfile=$1
  TMPFILE1=`mktemp DUSTOUT.XXXXXX` || exit 1
  time ${DUSTMASKER} -in ${inputfile} -linker 5 -outfmt fasta -level 15 -window 32 -out ${TMPFILE1}
  TMPFILE2=`mktemp ENCSEQ_DUST_OUT.XXXXXX` || exit 1
  time $2 bin/gt encseq encode -dust -dustecho -dustlink 5 -dustthreshold 1.5 -dustwindow 32 ${inputfile} > ${TMPFILE2}
  diff -I '^[>#]' ${TMPFILE1} ${TMPFILE2}
  rm -f ${TMPFILE2} ${TMPFILE1}
}

dust_compare_all()
{
  TMPFILE=`mktemp SEQ.XXXXXX` || exit 1
  for filename in `${GTDIR}/scripts/findfasta.rb --no-fastq --no-gzip --excludelist RandomN.fna`
  do
    echo ">$filename" > ${TMPFILE}
    egrep -v '^>' ${filename} | sed -e 's/[^ACGTacgt]//g' | tr acgt ACGT >> ${TMPFILE}
    dust_compare ${TMPFILE} $1
  done
  rm -f ${TMPFILE}
}

if test $# -eq 1
then
  if $1 == "valgrind"
  then
    dust_compare_all "valgrind"
  else
    dust_compare $1
  fi
else
  if test $# -eq 0
  then
    dust_compare_all
  else
    if test $# -eq 2
    then
      dust_compare $2 $1
    else
      echo "Usage: $0 [valgrind] [file]"
      exit 1
    fi
  fi
fi
