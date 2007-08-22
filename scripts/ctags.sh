#!/bin/sh
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
# See LICENSE file or http://genometools.org/license.html for license details.
#

# make a new tags file
ctags src/*.[ch]            \
      src/libgtcore/*.[ch]  \
      src/libgtext/*.[ch]   \
      src/libgtmatch/*.[ch] \
      src/libgtview/*.[ch]  \
