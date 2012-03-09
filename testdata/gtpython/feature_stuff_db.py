#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2008-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008-2012 Center for Bioinformatics, University of Hamburg
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


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: " + (sys.argv)[0] + " GFF3_db\n")
        sys.stderr.write('Test the FeatureIndex bindings on GFF3_db.\n')
        sys.exit(1)

  # try to instantiate an interface -> should fail

    try:
        feature_index = FeatureIndex()
    except NotImplementedError:
        pass
    else:
        raise

    rdb = RDBSqlite(sys.argv[1])
    adbs = AnnoDBGFFLike()

    feature_index = adbs.get_feature_index(rdb)
    seqid = feature_index.get_first_seqid()
    features = feature_index.get_features_for_seqid(seqid)
    if not features:
        raise TestFailedError

    rng = feature_index.get_range_for_seqid(seqid)
    features_r = feature_index.get_features_for_range(rng.start, rng.end,
            seqid)

    if not len(features) == len(features_r):
        raise TestFailedError

    gff3_visitor = GFF3Visitor()

    for feature in features:
        feature.accept(gff3_visitor)
        if not isinstance(feature, FeatureNode):
            raise TestFailedError
