#!/bin/sh

if test -f lib/libgenometools.a
then
  extractsyms.sh lib/libgenometools.a lib/libgtunstable.a | sort -u > GTNAMES
  common=`comm -1 -2 GTNAMES ${WORKVSTREE}/src/VMNAMES | wc -l`
  if test $common -gt 0
  then
    echo "========= the folllowing name clashes have been identified ======"
    comm -1 -2 GTNAMES ${WORKVSTREE}/src/lib/VMNAMES
  else
    echo "========= no name clashes have been identified ======"
  fi
  rm -f GTNAMES
else
  echo "========= no genometools lib available ======"
  exit 1
fi
