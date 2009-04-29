#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008 Center for Bioinformatics, University of Hamburg
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

from gt.core import *
from gt.extended import *
from gt.annotationsketch import *
import sys
import re

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: " + (sys.argv)[0] + " GFF3_file\n")
        sys.stderr.write("Show sequence ids contained in GFF3 annotation file.")
        sys.exit(1)

    in_stream = GFF3InStream((sys.argv)[1])
    feature_index = FeatureIndexMemory()
    feature_stream = FeatureStream(in_stream, feature_index)
    gn = feature_stream.next_tree()

  # fill feature index

    while gn:
        gn = feature_stream.next_tree()

    seqids = feature_index.get_seqids()
    for seqid in seqids:
        print seqid
