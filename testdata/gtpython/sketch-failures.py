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

from ctypes import ArgumentError
from gt.core import *
from gt.extended import *
from gt.annotationsketch import *
import sys
import re


class TestFailedError(Exception):

    pass


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: " + (sys.argv)[0] + " GFF3_file\n")
        sys.stderr.write("Test failures in GFF3 annotation file.")
        sys.exit(1)

    feature_index = FeatureIndexMemory()
    feature_index.add_gff3file((sys.argv)[1])

    seqid = feature_index.get_first_seqid()
    range = feature_index.get_range_for_seqid(seqid)

    style = Style()
    style.set_num("format", "margins", 40)

  # check error reporting in the Diagram class

    try:
        diagram = Diagram.from_index(feature_index, "nonexist", range,
                style)
    except GTError, strerr:
        if -1 == str(strerr).find("feature index does not contain "):
            raise TestFailedError
    else:
        raise TestFailedError
    try:
        diagram = Diagram.from_index(feature_index, seqid, range, 42)
    except ArgumentError, strerr:
        if -1 == str(strerr).find("must be a Style"):
            raise TestFailedError
    else:
        raise TestFailedError
    try:
        diagram = Diagram.from_index(feature_index, seqid, 42, style)
    except AttributeError, strerr:
        if -1 == str(strerr).find("object has no attribute 'start'"):
            raise TestFailedError
    else:
        raise TestFailedError
    try:
        diagram = Diagram.from_index(42, seqid, range, style)
    except ArgumentError, strerr:
        if -1 == str(strerr).find("must be a FeatureIndex"):
            raise TestFailedError
    else:
        raise TestFailedError
    diagram = Diagram.from_index(feature_index, seqid, range, style)

  # check error reporting in the Layout class

    try:
        layout = Layout(diagram, 70, style)
    except GTError, strerr:
        if -1 == str(strerr).find("layout width must at least be twice"):
            raise TestFailedError
    else:
        raise TestFailedError
    layout = Layout(diagram, 700, style)
    height = layout.get_height()

  # check error reporting in the CanvasCairoFile class

    try:
        canvas = CanvasCairoFile(12, 700, height, None)
    except ArgumentError, strerr:
        if -1 == str(strerr).find("must be a Style"):
            raise TestFailedError
    else:
        raise TestFailedError
    canvas = CanvasCairoFile(style, 700, height, None)
