#!/bin/sh
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

SKTOOLS='src/tools/gt_suffixerator.c
        src/tools/gt_sfxmap.c
        src/tools/gt_trieins.c
        src/tools/gt_mergeesa.c'

for filename in `ls ${SKTOOLS} src/libgtmatch/*.c`
do
  headfile="obj/src/libgtmatch/`basename ${filename} .c`"
  rm -f ${headfile}.o ${headfile}.d ${headfile}.splint
done

rm -rf indexdir
rm -f lib/libgtmatch.a
rm -f testsuite/TMP.[a-zA-Z0-9]*
