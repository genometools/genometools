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

from comment_node import *
from feature_node import *
from custom_stream import *
from custom_stream_example import *
from custom_visitor import *
from genome_node import *
from genome_stream import *
from gff3_in_stream import *
from gff3_out_stream import *
from gff3_visitor import *
from node_visitor import *
from region_node import *
from sequence_node import *
from add_introns_stream import *
from inter_feature_stream import *
from dup_feature_stream import *
from merge_feature_stream import *
from rdb import *
from anno_db import *

CommentNode.register(gtlib)
GenomeNode.register(gtlib)
RegionNode.register(gtlib)
SequenceNode.register(gtlib)
MetaNode.register(gtlib)
EOFNode.register(gtlib)
CustomStream.register(gtlib)
FeatureNode.register(gtlib)
FeatureNodeIterator.register(gtlib)
GFF3InStream.register(gtlib)
AddIntronsStream.register(gtlib)
InterFeatureStream.register(gtlib)
DuplicateFeatureStream.register(gtlib)
MergeFeatureStream.register(gtlib)
RDBSqlite.register(gtlib)
AnnoDBSchema.register(gtlib)
