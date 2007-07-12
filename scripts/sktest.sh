#!/bin/sh
#
# Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

# Current Problems: MKVTREESMAPDIR  must be set
# Current Problems: for icc compiled version add 
# /usr/local/zbh/intel/cc/9.1/lib

set -e -x

outoptions="-tis -lcp -suf -bwt"

# the make call normally used for development
cd testsuite
for options in `alloutputoptions.rb`
do
  cmpdbfile.sh ${options} -pl 1 -db ${ATK}
done
checkmapped.sh -db ${ATK} ${AT} ${GRUMBACH}/*.fna -parts 3 -pl 8 
checkmapped.sh -parts 1 -pl 3 -db ${SWK} ${SW}
checkmapped.sh -db ${SWK} ${SW} -parts 3 -pl 3 
checkmapped.sh -parts 2 -pl 7 -smap TransDNA -db ${AT}
checkmapped.sh -db ${SWK} -parts 1 -pl 3 -smap TransProt11 
cmpdbfile.sh ${outoptions} -pl 1 -db ../testdata/Random-Small.fna
cmpdbfile.sh ${outoptions} -pl 1 -db ../testdata/Random.fna
cmpdbfile.sh ${outoptions} -pl 1 -db ../testdata/Atinsert.fna ../testdata/Random.fna
cmpdbfile.sh ${outoptions} -pl 1 -db ../testdata/TTT-small.fna
cmpdbfile.sh ${outoptions} -pl 1 -db ${ATK} ${AT} ${GRUMBACH}/*.fna
testsuite.rb -keywords gt_suffixerator
cd ..
