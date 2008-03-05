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

# skproto-all.sh

makecompilerflags()
{
  printf "all:\n\t\${MAKE} curses=no"
  if test $3 -eq 1
  then
    printf " CC='ccache icc'"
  else
    printf " CC='ccache gcc'"
  fi
  if test $1 -eq 64
  then
    printf " 64bit=yes"
  fi
  printf " CFLAGS='-O3 -m$1"
  # printf " -DINLINEDENCSEQ"
  # printf " -DINLINEDSequentialsuffixarrayreader"
  if test $2 -eq 1
  then
    printf " -DBIGSEQPOS"
  fi
  if test $3 -eq 1
  then
    printf " -wd1418,869,981,1338"
  fi
  printf "'"
  printf " LDFLAGS='-m$1'"
  if test $3 -eq 1
  then
    printf " LD='icc' CXX='icc'"
  fi
  # printf " -DWITHTRIEIDENT"
  printf "\n"
}

# printf "NOASSERT='assert=no'"
if test -r LocalMakefile
then
  mv LocalMakefile LocalMakefile.previous
fi

if test $icc -eq 1
then
  makecompilerflags 32 $big $icc > LocalMakefile
else
  if test $do64 -eq 1
  then
    makecompilerflags 64 0 $icc > LocalMakefile
  else
    makecompilerflags 32 $big $icc > LocalMakefile
  fi
fi

if test -r LocalMakefile.previous
then
  cmp -s LocalMakefile LocalMakefile.previous
  if test $? -eq 1
  then
    echo "Current and previous LocalMakefile files differ: first remove them"
    exit 1
  fi
fi
make -f LocalMakefile

# echo ${MAKEFLAGS}
# make ${MAKEFLAGS} bin/skproto
# make ${MAKEFLAGS} $*

#make splint-gtmatch
