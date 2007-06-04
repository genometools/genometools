#!/bin/sh
#
# Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

set -e -x

# the make call normally used for development
cd testsuite
checkmapped.sh -parts 3 -pl 8 ${ATK} ${AT} ${GRUMBACH}/*.fna
checkmapped.sh -parts 1 -pl 3 ${SWK} ${SW}
checkmapped.sh -parts 3 -pl 3 ${SWK} ${SW}
checkmapped.sh -parts 2 -pl 7 -smap TransDNA ${AT}
checkmapped.sh -parts 1 -pl 3 -smap TransProt11 ${SWK}
cmpdbfile.sh ../testdata/Random-Small.fna
cmpdbfile.sh ../testdata/Random.fna
cmpdbfile.sh ../testdata/Atinsert.fna ../testdata/Random.fna
cmpdbfile.sh ../testdata/TTT-small.fna
cmpdbfile.sh ${ATK} ${AT} ${GRUMBACH}/*.fna
testsuite.rb -keywords gt_suffixerator
cd ..
