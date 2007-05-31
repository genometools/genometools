#!/bin/sh
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

# the make which links cairo
gmake CFLAGS=-I/usr/local/include/cairo LDFLAGS='-L/usr/local/lib -L/usr/X11R6/lib' LDLIBS='-lm -lz -lcairo' $*
