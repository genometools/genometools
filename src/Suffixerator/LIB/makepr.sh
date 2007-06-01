#!/bin/sh

for filename in `ls *.c`
do
  skproto.x ${filename} > `basename ${filename} .c`.pr
done
