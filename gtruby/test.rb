#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
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

GT::allocators_init()
GT::allocators_reg_atexit_func()

error = GT::Error.new()

in_stream = GT::GFF3InStream.new("../testdata/encode_known_genes_Mar07.gff3")
#in_stream = GT::GFF3InStream.new("../testdata/standard_gene_as_tree.gff3")
#out_stream = GT::GFF3OutStream.new(in_stream)
gff3_visitor = GT::GFF3Visitor.new()
gn = in_stream.next_tree()
gn.accept(gff3_visitor)
while gn do
  gn = in_stream.next_tree()
  if gn then gn.accept(gff3_visitor) end
end



=begin
feature_index = GT::FeatureIndex.new()
feature_stream = GT::FeatureStream.new(in_stream, feature_index)
while gn do
  gn = feature_stream.next_tree()
end
=end
