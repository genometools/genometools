#!/usr/bin/env ruby
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

require 'gtruby'

if ARGV.size != 3 then
  STDERR.puts "Usage: #{$0} style_file PNG_file GFF3_file"
  STDERR.puts "Create PNG representation of GFF3 annotation file."
  exit(1)
end

(stylefile, pngfile, gff3file) = ARGV

# load style file
style = GT::Style.new()
style.load_file(stylefile)

# create feature index
feature_index = GT::FeatureIndexMemory.new()

# add GFF3 file to index
feature_index.add_gff3file(gff3file)

# create diagram for first sequence ID in feature index
seqid = feature_index.get_first_seqid()
range = feature_index.get_range_for_seqid(seqid)
diagram = GT::Diagram.from_index(feature_index, seqid, range, style)

# create layout for given width
layout = GT::Layout.new(diagram, 800, style)

# create canvas with given width and computed height
canvas = GT::CanvasCairoFile.new(style, 800, layout.get_height, nil)

# sketch layout on canvas
layout.sketch(canvas)

# write canvas to file
canvas.to_file(pngfile)
