#!/bin/sh
#
# Copyright (c) 2004-2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2004-2006 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

# clean (i.e., remove!) all files comprising the index constructed by mkvtree
# and all .polya and .polya.info files constructed by gth
# furthermore, all bioseq files (.gt_bsi .gt_bsr) are removed

rm -f *.ssp *.llv *.skp *.al[12] *.des *.prj *.lcp *.suf
rm -f *.tis *.ois *.bwt *.bck *.sds *.iso *.sti *.sti1 *.cld
rm -f *.cld1 *.crf *.cfr *.lsf
rm -f *.polya *.polya.info
rm -f *.gt_bsi *.gt_bsr
