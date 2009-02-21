#!/bin/sh

VSTREEBIN=/projects/vstree/bin/i686-pc-linux-gnu

if test ! -d ${VSTREEBIN}
then
  VSTREEBIN=/Users/stefan/vstree/bin/i686-apple-darwin
fi

set  -e -x

if test $# -ne 4
then
  echo "Usage: $0 <number of lines> <# refseq> <# query seq> <threshold>"
  exit 1
fi

# modify the following if you need larger files.
# but only use an extractfile that does not contain wildcard characters,
# i.e. N etc.

extractfile=testdata/U89959_genomic.fas
lines=$1
seqref=$2
seqquery=$3
threshold=$4

createseq()
{
  TMPFILE=`mktemp gt-test.XXXXXX` || exit 1

  echo ">$3 -n $1" > ${TMPFILE}.extract
  grep -v '^>' ${extractfile} | $3 -n $1 >> ${TMPFILE}.extract
  gt mutate -rate 10 ${TMPFILE}.extract | grep -v '^#' > ${TMPFILE}.mut
  ${VSTREEBIN}/mkvtree -db ${TMPFILE}.mut -indexname ${TMPFILE} -ois -tis -dna
  ${VSTREEBIN}/vsubseqselect -minlength 100 -maxlength 200 -snum $2 ${TMPFILE} |\
                  grep -v '^#' |
                  gawk '/.*/ {print ">\n" $1}'
  rm -f ${TMPFILE}.*
}

createseq $lines $seqref head > gt-test-refer
createseq $lines $seqquery tail > gt-test-query

scripts/failure.sh gt-test-refer gt-test-query ${threshold}
