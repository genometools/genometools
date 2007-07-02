#!/bin/bash
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

for filename in `ls src/*.c src/libgtmatch/*.c`
do
  headfile="obj/`basename ${filename} .c`"
  rm -f ${headfile}.o ${headfile}.d ${headfile}.splint
done

rm -f lib/libgtmatch.a
rm -f testsuite/TMP.[A-Z0-9]*
