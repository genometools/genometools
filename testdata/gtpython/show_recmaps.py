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
        sys.stderr.write("Show RecMaps of GFF3 annotation file.")
        sys.exit(1)

    in_stream = GFF3InStream((sys.argv)[1])
    feature_index = FeatureIndexMemory()
    feature_stream = FeatureStream(in_stream, feature_index)
    gn = feature_stream.next_tree()

  # fill feature index

    while gn:
        gn = feature_stream.next_tree()

    seqid = feature_index.get_first_seqid()
    range = feature_index.get_range_for_seqid(seqid)

    style = Style()
    diagram = Diagram.from_index(feature_index, seqid, range, style)
    layout = Layout(diagram, 700, style)
    height = layout.get_height()
    image_info = ImageInfo()
    canvas = CanvasCairoFile(style, 800, height, image_info)
    layout.sketch(canvas)
    canvas.to_file("foo1.png")

    for (x1, y1, x2, y2, gn) in image_info.each_hotspot():
        print("x1=%d, y1=%d, x2=%d, y2=%d, gn.type=%s" % (x1, y1, x2, y2,
                gn.get_type()))
