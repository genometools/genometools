#!/bin/sh

# author: Stefan Kurtz, Nov. 2015

set -e -x

if test $# -ne 2
then
  echo "Usage: $0 <minlen> <inputfile>"
  exit 1
fi

minlen=$1
inputfile=$2
TARGET=`basename $inputfile`
minidentity=70
seedlength=14
maxfreq=21

# note that the following script randomly replaces wildcards by
# characters over the base alphabet. So for the same sequence and parameters
# the set of matches often varies over different runs.
scripts/convert2myersformat.rb $inputfile > ${TARGET}.fasta
if test $PACKAGES != ""
then
  MYERSPROG="$PACKAGES"/myers
else
  MYERSPROG="../myers"
fi

rm -f ${TARGET}.db
${MYERSPROG}/DAZZ_DB/fasta2DB ${TARGET}.db ${TARGET}.fasta
${MYERSPROG}/DALIGNER/daligner -t${maxfreq} -I -A -Y -e0.${minidentity} \
                   -k${seedlength} -l${minlen} \
                   ${TARGET}.db ${TARGET}.db > ${TARGET}-da.matches
rm -f ${TARGET}.db
rm -f ${TARGET}.${TARGET}*.las .${TARGET}.idx .${TARGET}.bps
bin/gt encseq encode -sds no -md5 no -des no -indexname ${TARGET} ${TARGET}.fasta
bin/gt seed_extend -ii ${TARGET} -t ${maxfreq} -l ${minlen} \
                    -minidentity ${minidentity} -seed-display -v \
                    -overlappingseeds -bias-parameters -history 60 > ${TARGET}-se.matches
scripts/matched-seqpairs.rb ${TARGET}-da.matches ${TARGET}-se.matches
