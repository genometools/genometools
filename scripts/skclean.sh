#!/bin/bash
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

localpath=src/libgtmatch

for filename in `ls ${localpath}/*.c`
do
  objfile="obj/`basename ${filename} .c`.o"
  rm -f $objfile
done

rm -f lib/libgtmatch.a
