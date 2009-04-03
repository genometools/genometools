#!/bin/sh

if test -f lib/libgenometools.a
then
  extractsyms.sh lib/libgenometools.a lib/libgtunstable.a | sort -u > GTNAMES
  echo "========= the folllowing name clashes have been identified ======"
  comm -1 -2 GTNAMES ${WORKVSTREE}/src/lib/VMNAMES
  rm -f GTNAMES
else
  echo "========= no genometools lib available ======"
  exit 1
fi
