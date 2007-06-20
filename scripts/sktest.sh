#!/bin/sh
#
# Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

set -e -x

# the make call normally used for development
cd testsuite
checkmapped.sh -db ${ATK} ${AT} ${GRUMBACH}/*.fna -parts 3 -pl 8 
checkmapped.sh -parts 1 -pl 3 -db ${SWK} ${SW}
checkmapped.sh -db ${SWK} ${SW} -parts 3 -pl 3 
checkmapped.sh -parts 2 -pl 7 -smap TransDNA -db ${AT}
checkmapped.sh -db ${SWK} -parts 1 -pl 3 -smap TransProt11 
cmpdbfile.sh ../testdata/Random-Small.fna
cmpdbfile.sh ../testdata/Random.fna
cmpdbfile.sh ../testdata/Atinsert.fna ../testdata/Random.fna
cmpdbfile.sh ../testdata/TTT-small.fna
cmpdbfile.sh ${ATK} ${AT} ${GRUMBACH}/*.fna
testsuite.rb -keywords gt_suffixerator
cd ..
