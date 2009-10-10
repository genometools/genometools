#!/bin/sh
#
# Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
# Copyright (c) 2009 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

set -e -x

outoptions="-tis -lcp -suf -bwt"

cd testsuite
num=2
while test ${num} -lt 10 
do
  ../scripts/iterrunmerge.sh ${num}
  num=`expr ${num} + 1`
done

#depends on mkvtree.x

../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/Random-Small.fna
../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/Random.fna
../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/Atinsert.fna ../testdata/Random.fna
../scripts/cmpdbfile.sh ${outoptions} -pl -db ../testdata/TTT-small.fna
ALLOUTPUTOPTS="../scripts/alloutputoptions.rb"
if test ! -f ${ALLOUTPUTOPTS}
then
  echo "cannot find ${ALLOUTPUTOPTS}"
  exit 1
fi

if test ! "X${GTTESTDATA}" = "X"
then
  AT=../testdata/at1MB
  ATK=${GTTESTDATA}/Iowa/at100K1
  GRUMBACH=${GTTESTDATA}/DNA-mix/Grumbach.fna
  for options in `${ALLOUTPUTOPTS}`
  do
    ../scripts/cmpdbfile.sh ${options} -pl -db ${ATK}
  done
  ../scripts/cmpdbfile.sh ${outoptions} -pl -db ${ATK} ${AT} ${GRUMBACH}/*.fna
fi

if test ! "X${GTTESTDATA}" = "X"
then
  AT=../testdata/at1MB
  U8=../testdata/U89959_genomic.fas
  ../scripts/rununique.sh ../bin/gt 10 20 ${U8} ${AT}
fi

cd ..
