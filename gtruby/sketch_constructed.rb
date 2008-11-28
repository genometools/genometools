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

if ARGV.size != 2 then
  STDERR.puts "Usage: #{$0} style_file PNG_file"
  exit(1)
end

seqid = "chromosome_21"

# construct a gene on the forward strand with two exons
gene   = GT::FeatureNode.create(seqid, "gene", 100, 900, "+")
exon   = GT::FeatureNode.create(seqid, "exon", 100, 200, "+")
gene.add_child(exon)
intron = GT::FeatureNode.create(seqid, "intron", 201, 799, "+")
gene.add_child(intron)
exon   = GT::FeatureNode.create(seqid, "exon", 800, 900, "+")
gene.add_child(exon)

# construct a single-exon gene on the reverse strand
# (within the intron of the forward strand gene)
reverse_gene = GT::FeatureNode.create(seqid, "gene", 400, 600, "-")
reverse_exon = GT::FeatureNode.create(seqid, "exon", 400, 600, "-")
reverse_gene.add_child(reverse_exon)

pngfile = ARGV[1]

style = GT::Style.new()
style.load_file(ARGV[0])

rng = GT::Range.malloc()
rng.start = 1
rng.end   = 1000

diagram = GT::Diagram.from_array([gene, reverse_gene], rng, style)

layout = GT::Layout.new(diagram, 600, style)
canvas = GT::CanvasCairoFile.new(style, 600, layout.get_height, nil)
layout.sketch(canvas)

canvas.to_file(pngfile)
