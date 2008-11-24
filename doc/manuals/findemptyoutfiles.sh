#!/bin/sh

for filename in `find output/ -name '*.out'`
do
  filesize=`cat ${filename} | wc -l`
  if test $filesize -eq 1
  then
	echo "Warning: ${filename} only contains one line"
        exit 1
  fi
done
