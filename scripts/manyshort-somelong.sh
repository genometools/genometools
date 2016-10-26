#!/bin/sh

indexname=manyshort-onelong.fna

if test $# -ne 2
then
  echo "Usage: $0 <directory> <numofsequences>"
  exit 1
fi

directory=$1
numofsequences=$2

for inputfile in `ls ${directory}/*.fna`
do
  scripts/cutsequences.rb ${inputfile} $numofsequences 50 50
  cat ${inputfile}
done
