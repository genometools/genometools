#!/bin/sh
if test $# -ne 3
then
  echo "Usage: $0 <readsize> <coverage> <inputfile>"
  exit 1
fi

readsize=$1
coverage=$2
inputfile=$3

if test -f bin/gt
then
  GTBIN=bin/gt
else
  GTBIN=${GTINSTALL}/bin/gt
fi

${GTBIN} suffixerator -v -des no -sds no -tis -ssp -dna -db ${inputfile} -indexname genome-idx
${GTBIN} simreads -coverage ${coverage} -len ${readsize} -gzip -force \
            -o genome-idx-${readsize}-${coverage}-reads.fna.gz genome-idx
rm -f genome-idx.*
