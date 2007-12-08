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

# testing some Ruby bindings for FeatureIndex

require 'gtruby'

if ARGV.size != 1 then
  STDERR.puts "Usage: #{$0} GFF3_file"
  STDERR.puts "Show sequence ids contained in GFF3 annotation file."
  exit(1)
end


gff3file = ARGV[0]

in_stream = GT::GFF3InStream.new(gff3file)
feature_index = GT::FeatureIndex.new()
feature_stream = GT::FeatureStream.new(in_stream, feature_index)
gn = feature_stream.next_tree()
# fill feature index
while (gn) do
  gn = feature_stream.next_tree()
end

seqids = feature_index.get_seqids()
seqids.each { |seqid| puts seqid }
