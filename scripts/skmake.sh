#!/bin/bash
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

localpath=src/libgtmatch
usage="Usage: $0 [icc,64]"

if test $# -eq 0
then
  icc=0
  do64=0
else
  if test $# -eq 1
  then
    case $1 in
     "icc") icc=1;;
     "64")  do64=1;;
     *) echo $usage
        exit 1;;
    esac
  else
    echo $usage
    exit 1
  fi
fi

for filename in `ls ${localpath}/*.c`
do
  prfile="${localpath}/`basename ${filename} .c`.pr"
  if test ${filename} -nt ${prfile}
  then
    skproto.x ${filename} > ${prfile}
  fi
done

if test $icc -eq 1
then
  make CC='icc' CFLAGS='-O3 -wd1418,869,981' LD='icc' CXX='icc'
else
  if test $do64 -eq 1
  then
    make CC='ccache gcc' CFLAGS='-O3 -m64' LDFLAGS='-m64'
  else
    make CC='ccache gcc' CFLAGS='-O3 -m32' LDFLAGS='-m32'
  fi
fi
