#!/bin/sh

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
rm -f ${TARGET}.db *.las
bin/gt encseq encode ${TARGET}.fasta
bin/gt seed_extend -ii ${TARGET}.fasta -maxfreq ${maxfreq} -l ${minlen} \
                    -minidentity ${minidentity} -seed-display -v \
                    -overlappingseeds > ${TARGET}-sa.matches
