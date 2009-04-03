#!/bin/sh

#set -e -x

if test $# -eq 0
then
  echo "Usage: $0 <archive1> <archive2> ..."
  exit 1
fi

for archive in $*
do
  ar x $archive
  for filename in `ls *.o`
  do
    nm ${filename} | gawk '/^[0-9]/ {if ($2 == "t" || $2 == "T") print $3}'
    rm ${filename}
  done
done
