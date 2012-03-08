#
# Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2009 Center for Bioinformatics, University of Hamburg
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
  STDERR.puts "Test the custom stream bindings on GFF3 file."
  exit(1)
end

module GT
  class PrinterStream < CustomStream
    def initialize(instream, types)
      super()
      @types = types
      @instream = instream
      @visitor = PrinterVisitor.new(@types)
    end

    def next
      node = @instream.next_tree
      node.accept(@visitor) if node
      node
    end
  end

  class PrinterVisitor < CustomVisitor
    def initialize(types)
      super()
      @types = types
    end

    def visit_feature_node(fn)
      fn.traverse_dfs do |node|
        @types.push(node.get_type)
        @types.uniq!
        rng = node.get_range
        puts "#{node.get_source} #{node.get_type} #{rng.begin} #{rng.end}"
        attribs = []
        node.each_attribute do |k,v|
           attribs.push("#{k}=#{v}")
        end
        puts attribs.sort.join(',')
      end
      0
    end
  end

  class ErrorStream < CustomStream
    def initialize
      super
    end

    def next
      raise "foo"
    end
  end

  class IncompleteStream < CustomStream
    def initialize
      super
    end
  end
end

gff3file = ARGV[0]

types = []
is = GT::GFF3InStream.new(gff3file)
ps = GT::PrinterStream.new(is, types)
while ps.next_tree do
end

raise unless types.include?("gene")
raise unless types.include?("mRNA")
raise unless types.include?("TF_binding_site")
raise unless types.include?("exon")
raise unless types.include?("CDS")

es = GT::ErrorStream.new()
begin
  while es.next_tree do
  end
rescue GT::GTError
  # exception expected
else
  raise
end

begin
  is = GT::IncompleteStream.new
rescue GT::GTError
  # exception expected
else
  raise
end
