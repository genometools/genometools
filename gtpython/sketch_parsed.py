#!/usr/bin/env python
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

from gt.annotationsketch import *
from gt.core.gtrange import Range
import sys

if __name__ == "__main__":
  if len(sys.argv) != 4:
    sys.stderr.write("Usage: " + sys.argv[0] + " Style_file PNG_file GFF3_file\n")
    sys.stderr.write("Create PNG representation of GFF3 annotation file.")
    sys.exit(1)

  pngfile = sys.argv[2]

  # load style file
  style = Style()
  style.load_file(sys.argv[1])

  # create feature index
  feature_index = FeatureIndexMemory()

  # add GFF3 file to index
  feature_index.add_gff3file(sys.argv[3])

  # create diagram for first sequence ID in feature index
  seqid = feature_index.get_first_seqid()
  range = feature_index.get_range_for_seqid(seqid)
  diagram = Diagram.from_index(feature_index, seqid, range, style)

  # create layout
  layout = Layout(diagram, 600, style)
  height = layout.get_height()

  # create canvas
  canvas = CanvasCairoFile(style, 600, height)

  # sketch layout on canvas
  layout.sketch(canvas)

  # write canvas to file
  canvas.to_file(pngfile)
