#!/bin/sh
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

# the make call normally used for development
make CC='gcc' CFLAGS='-O3 -m32' LDFLAGS='-m32'
