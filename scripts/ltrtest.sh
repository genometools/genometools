#!/bin/sh
#
# Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

set -e -x

# the make call normally used for development
cd testsuite
env -i ./testsuite.rb -keywords gt_ltr
env -i ./testsuite.rb -keywords 'gt_ltrharvest' -gttestdata ${GTTESTDATA}
cd ..
# optional -memcheck   (run valgrind)
#          -select 253 (run testcase 253)
#cd ../
