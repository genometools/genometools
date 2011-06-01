#!/bin/bash

TEMPLATE=src/match/esa-bottomup

for suffix in maxpairs shulen
do
  file=${TEMPLATE}-${suffix}.inc
  rm -f $file
  scripts/appendsuffix.rb ${TEMPLATE}.gen ${suffix} BUinfo BUstate > $file
done
