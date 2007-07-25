#!/bin/sh
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

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

skproto-all.sh

if test $big -eq 1
then
  bignum=-DBIGSEQPOS
else
  bignum=
fi

#-DWITHTRIEIDENT
FASTDEF="-Duint_fast32_t=uint32_t  -Duint_fast64_t=uint64_t"

if test $icc -eq 1
then
  make assert=no CC='ccache icc' CFLAGS='-O3 ${FASTDEF} ${bignum} -wd1418,869,981,1338' LD='icc' CXX='icc' $*
else
  if test $do64 -eq 1
  then
    make assert=no CC="ccache gcc" CFLAGS="-O3 -m64 ${FASTDEF}" LDFLAGS="-m64" $*
  else
    make assert=no CC="ccache gcc" CFLAGS="-O3 -m32 ${bignum} ${FASTDEF}" LDFLAGS="-m32" $*
  fi
fi

#make splint-gtmatch
