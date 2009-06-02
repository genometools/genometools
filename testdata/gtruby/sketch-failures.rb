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

if ARGV.size != 1 then
  STDERR.puts "Usage: #{$0} GFF3_file"
  exit(1)
end

class TestFailedError < Exception
end

feature_index = GT::FeatureIndexMemory.new()
feature_index.add_gff3file(ARGV[0])

seqid = feature_index.get_first_seqid()
range = feature_index.get_range_for_seqid(seqid)

style = GT::Style.new()
style.set_num("format","margins",100)

begin
  diagram = GT::Diagram.from_index(feature_index, "nonexist", range, style)
rescue GT::GTError => msg
    raise if !/feature index does not contain the given sequence id/.match(msg)
else
  raise TestFailedError
end
begin
  diagram = GT::Diagram.from_index(nil, seqid, range, style)
rescue
  # exception expected
else
  raise TestFailedError
end
begin
  diagram = GT::Diagram.from_index(feature_index, seqid, nil, style)
rescue
  # exception expected
else
  raise TestFailedError
end
begin
  diagram = GT::Diagram.from_index(feature_index, seqid, range, "Dd")
rescue
  # exception expected
else
  raise TestFailedError
end
diagram = GT::Diagram.from_index(feature_index, seqid, range, style)

begin
  layout = GT::Layout.new(diagram, 70, style)
rescue GT::GTError => msg
  raise if !/layout width must at least be twice the x-margin/.match(msg)
else
  raise TestFailedError
end
layout = GT::Layout.new(diagram, 700, style)

canvas = GT::CanvasCairoFile.new(style, 700, layout.get_height, nil)
