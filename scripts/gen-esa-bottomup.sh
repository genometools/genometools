#!/bin/bash

TEMPLATE=src/match/esa-bottomup

for suffix in maxpairs shulen spmsk
do
  file=${TEMPLATE}-${suffix}.inc
  rm -f $file
  scripts/appendsuffix.rb ${TEMPLATE}.gen ${suffix} GtBUinfo GtBUstate \
                          GtArrayGtBUItvinfo > $file
done
