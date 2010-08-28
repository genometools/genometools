#!/bin/sh

splitmultifasta.rb tmp 10 ../testdata/Atinsert.fna
for query in `ls tmp-[0-9]`
do
  for refer in `ls tmp-[0-9]`
  do
    shulength.sh ${query} ${refer}
    if test $? -ne 0
    then
      echo "failure: shulength.sh ${query} ${refer}"
      exit 1
    fi
  done
done
