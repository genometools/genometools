#!/bin/sh

if test $# -eq 0
then
  echo "Usage: missing arguments"
  exit 1
fi

#TMPFILE1=`mktemp /tmp/maketags.XXXXXX` || exit 1
#TMPFILE2=`mktemp /tmp/maketags.XXXXXX` || exit 1
 
# grep -v '^>' | tr -d '\n' | tr 'BDEFHIJKLMNOPQRSUVWXYZ:e'\
mkvtree.x -db $* -tis -ois -dna -indexname selectidx
vsubseqselect -minlength 32 -maxlength 32 -snum 10000 selectidx |\
                             grep '^#' -v |
                             sort |\
                             awk '{print ">\n" $1}'
# zcat $* | grep -v '^>' | tr -d '\n' > ${TMPFILE1}
# onePE.rand.x ${TMPFILE1} 1000 32 | sort | awk '{print ">\n" $1}' >  ${TMPFILE2}
# env -i GT_MEM_BOOKKEEPING=off bin/gt mutate -rate 1 ${TMPFILE2}
#rm -f ${TMPFILE1} ${TMPFILE2}
