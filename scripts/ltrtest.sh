#!/bin/sh
#
# Copyright (c) 2007 David Ellinghaus <dellinghaus@stud.zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

set -e -x

# the make call normally used for development
cd ../testsuite
./testsuite.rb -keywords gt_ltr
# optional -memcheck   (run valgrind)
#          -select 253 (run testcase 253)
#cd ../
