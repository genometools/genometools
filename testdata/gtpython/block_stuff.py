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


class TestFailedError(Exception):

    pass


def testfunc(bl):
    b = Block(bl)
    print "%s %s %s %d" % (b.get_type(), b.get_strand(), b.get_range(),
                           b.get_size())
    if not b.get_top_level_feature():
        raise TestFailedError
    return b.get_type()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: " + (sys.argv)[0] + " GFF3_file\n")
        sys.stderr.write("Test the TrackSelector callback and Block bindings.")
        sys.exit(1)

    genome_stream = GFF3InStream((sys.argv)[1])

  # instantiate index object

    feature_index = FeatureIndexMemory()
    genome_stream = FeatureStream(genome_stream, feature_index)

    feature = genome_stream.next_tree()
    while feature:
        feature = genome_stream.next_tree()

    seqid = feature_index.get_first_seqid()
    range = feature_index.get_range_for_seqid(seqid)

    style = Style()
    diagram = Diagram.from_index(feature_index, seqid, range, style)
    diagram.set_track_selector_func(testfunc)

    layout = Layout(diagram, 800, style)
    height = layout.get_height()

