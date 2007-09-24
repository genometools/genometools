#!/bin/sh
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
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

#skproto-all.sh

if test $big -eq 1
then
  bignum=-DBIGSEQPOS
else
  bignum=
fi

# NOASSERT='assert=no'
#-DWITHTRIEIDENT
# EXTRAFLAGS="-Duint_fast32_t=uint32_t  -Duint_fast64_t=uint64_t"
EXTRAFLAGS="-DINLINEDENCSEQ -DINLINEDSequentialsuffixarrayreader"
# EXTRAFLAGS="-DINLINEDENCSEQ"
COMMON="curses=no"

cd ../
if test $icc -eq 1
then
  make ${COMMON} CC='ccache icc' CFLAGS='-O3 ${EXTRAFLAGS} ${bignum} -wd1418,869,981,1338' LD='icc' CXX='icc' $*
else
  if test $do64 -eq 1
  then
    make ${COMMON} CC="ccache gcc" CFLAGS="-O3 -m64 ${EXTRAFLAGS}" LDFLAGS="-m64" $*
  else
    make ${COMMON} CC="ccache gcc" CFLAGS="-O3 -m32 ${bignum} ${EXTRAFLAGS}" LDFLAGS="-m32" $*
  fi
fi
cd ./scripts

#make splint-gtmatch
