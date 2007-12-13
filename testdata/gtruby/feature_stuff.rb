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

# testing the Ruby bindings for FeatureIndex and FeatureStream classes

require 'gtruby'

if ARGV.size != 1 then
  STDERR.puts "Usage: #{$0} GFF3_file"
  STDERR.puts "Test the FeatureIndex and FeatureStream bindings on GFF3 file."
  exit(1)
end


gff3file = ARGV[0]

# set up the feature stream
genome_stream = GT::GFF3InStream.new(gff3file)
feature_index = GT::FeatureIndex.new()
genome_stream = GT::FeatureStream.new(genome_stream, feature_index)

feature = genome_stream.next_tree()
while (feature) do
  feature = genome_stream.next_tree()
end

seqid = feature_index.get_first_seqid()
features = feature_index.get_features_for_seqid(seqid)
raise if not features
gff3_visitor = GT::GFF3Visitor.new()

features.each do |feature|
  feature.accept(gff3_visitor)
end
