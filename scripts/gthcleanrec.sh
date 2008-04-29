#!/bin/sh -e
#
# Copyright (c) 2004-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg
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

# descend recursively in all subdirectries of the current directory and
# clean (i.e., remove!) all files comprising the index constructed by mkvtree
# and all .polya and .polya.info files constructed by gth
# furthermore, all bioseq files (.gt_bsi .gt_bsr) are removed

SUFFIXES="ssp llv skp al1 al2 des prj lcp suf tis ois bwt bck sds iso sti sti1 cld cld1 crf cfr lsf polya polya.info gt_bsi gt_bsr gt_bsf"

for SUFFIX in $SUFFIXES
do
  find . -name "*.$SUFFIX" -print0 | xargs -0 rm -f
done
