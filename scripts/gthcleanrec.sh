#!/bin/sh -e
#
# Copyright (c) 2004-2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2004-2006 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

# descend recursively in all subdirectries of the current directory and
# clean (i.e., remove!) all files comprising the index constructed by mkvtree
# and all .polya and .polya.info files constructed by gth
# furthermore, all bioseq files (.gt_bsi .gt_bsr) are removed

SUFFIXES="ssp llv skp al1 al2 des prj lcp suf tis ois bwt bck sds iso sti sti1 cld cld1 crf cfr lsf polya polya.info gt_bsi gt_bsr"

for SUFFIX in $SUFFIXES
do
  find . -name "*.$SUFFIX" -print0 | xargs -0 rm -f
done
