#!/bin/bash
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

localpath=src/libgtmatch

for filename in `ls ${localpath}/*.c`
do
  prfile="${localpath}/`basename ${filename} .c`.pr"
  if test ${filename} -nt ${prfile}
  then
    skproto.x ${filename} > ${prfile}
  fi
done

# the make call normally used for development
make CC='ccache gcc' CFLAGS='-O3 -m32' LDFLAGS='-m32'
