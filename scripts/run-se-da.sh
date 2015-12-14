#!/bin/sh

# author: Stefan Kurtz, Nov. 2015

set -e -x

if test $# -eq 3
then
  minlen=$1
  minidentity=$2
  referencefile=$3
  queryfile=""
else
  if test $# -eq 4
  then
    minlen=$1
    minidentity=$2
    referencefile=$3
    queryfile=$4
  else
    echo "Usage: $0 <minlen> <minidentity> <referencefile> [queryfile]"
    exit 1
  fi
fi

reference=`basename $referencefile`
seedlength=14
maxfreq=21

# note that the following script randomly replaces wildcards by
# characters over the base alphabet. So for the same sequence and parameters
# the set of matches often varies over different runs.
scripts/convert2myersformat.rb $referencefile > ${reference}.fasta

if test $PACKAGES != ""
then
  MYERSPROG="$PACKAGES"/myers
else
  MYERSPROG="../myers"
fi

rm -f ${reference}.db
${MYERSPROG}/DAZZ_DB/fasta2DB ${reference}.db ${reference}.fasta
if test ${queryfile} != ""
then
  query=`basename $queryfile`
  scripts/convert2myersformat.rb ${queryfile} > ${query}.fasta
  ${MYERSPROG}/DAZZ_DB/fasta2DB ${query}.db ${query}.fasta
  outputprefix="${reference}-${query}"
else
  query=${reference}
  outputprefix="${reference}"
fi

${MYERSPROG}/DALIGNER/daligner -t${maxfreq} -I -A -Y -e0.${minidentity} \
                   -k${seedlength} -l${minlen} \
                   ${query}.db ${reference}.db > ${outputprefix}-da.matches
rm -f ${reference}.db ${query}.db .${reference}.idx .${reference}.bps
rm -f ${query}.${reference}*.las .${query}.idx .${query}.bps

bin/gt encseq encode -sds no -md5 no -des no -indexname ${reference} \
                                           ${reference}.fasta
if test ${queryfile} != ""
then
  bin/gt encseq encode -sds no -md5 no -des no -indexname ${query} \
                                           ${query}.fasta
  qiioption="-qii ${query}"
else
  qiioption=""
fi

bin/gt seed_extend -ii ${reference} ${qiioption} -t ${maxfreq} -l ${minlen} \
                    -seedlength 14 -minidentity ${minidentity} -seed-display \
                    -v -overlappingseeds -bias-parameters -no-reverse \
                    -history 60 > ${outputprefix}-se.matches
rm -f ${query}.ssp ${query}.fasta ${query}.esq
rm -f ${reference}.ssp ${reference}.fasta ${reference}.esq
scripts/matched-seqpairs.rb ${outputprefix}-da.matches ${outputprefix}-se.matches
