#!/usr/bin/env bash
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

localpath=src/libgtmatch
icc=0
do64=0
big=0

if test $# -ge 1
then
  case $1 in
   "-icc") icc=1
           shift;;
   "-m64")  do64=1
            shift;;
   "-big")  big=1
            shift;;
  esac
fi

for filename in `ls ${localpath}/*.c`
do
  prfile="${localpath}/`basename ${filename} .c`.pr"
  if test ${filename} -nt ${prfile}
  then
    bin/gt dev skproto ${filename} > ${prfile}
  fi
done

if test $big -eq 1
then
  bignum=-DBIGSEQPOS
else
  bignum=
fi

#-DWITHTRIEIDENT

if test $icc -eq 1
then
  make CC='ccache icc' CFLAGS='-O3 ${bignum} -wd1418,869,981,1338' LD='icc' CXX='icc' $*
else
  if test $do64 -eq 1
  then
    make CC="ccache gcc" CFLAGS="-O3 -m64" LDFLAGS="-m64" $*
  else
    make CC="ccache gcc" CFLAGS="-O3 -m32 ${bignum}" LDFLAGS="-m32" $*
  fi
fi

#make splint-gtmatch
