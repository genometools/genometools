#!/bin/bash

PRIMES=/projects/gi/primegen/primes
OUTFILE=`dirname $0`/../src/extended/md5set.prtab

echo "outfile: $OUTFILE"

function printprimes {
  echo "interval: [$1..$2]"
  echo "minimal distance: $3"
  $PRIMES $1 $2 $3 >> $4
  echo "[$1..$2] done"
}

mb1=1048576
mb4=4194304
mb64=67108864
gb1=1073741824
gb4=4294967296
gb64=68719476736
gb512=549755813888

rm $OUTFILE
touch $OUTFILE

#           min    max    step
printprimes $mb1   $gb4   $mb4   $OUTFILE
printprimes $gb4   $gb64  $mb64  $OUTFILE
printprimes $gb64  $gb512 $gb1   $OUTFILE

LARGEST_PRIME=`tail -n 1 $OUTFILE`
# remove trailing comma:
LARGEST_PRIME=${LARGEST_PRIME:0:${#LARGEST_PRIME}-1}

NOFPRIMES=`cat $OUTFILE | wc -l`

echo "#define GT_MD5SET_NOFPRIMES $NOFPRIMES" >> $OUTFILE
echo "#define GT_MD5SET_LARGEST_PRIME $LARGEST_PRIME" >> $OUTFILE
