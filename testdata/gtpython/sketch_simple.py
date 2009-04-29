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
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: " + (sys.argv)[0] +
                         " PNG_file GFF3_file\n")
        sys.stderr.write("Create PNG representation of GFF3 annotation file.")
        sys.exit(1)

    pngfile = (sys.argv)[1]
    gff3file = (sys.argv)[2]
    feature_index = FeatureIndexMemory()
    feature_index.add_gff3file(gff3file)

    seqid = feature_index.get_first_seqid()
    range = feature_index.get_range_for_seqid(seqid)

    style = Style()
    diagram = Diagram.from_index(feature_index, seqid, range, style)
    layout = Layout(diagram, 700, style)
    height = layout.get_height()
    canvas = CanvasCairoFile(style, 700, height)
    layout.sketch(canvas)
    canvas.to_file(pngfile)
