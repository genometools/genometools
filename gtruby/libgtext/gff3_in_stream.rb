#
# Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
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

require 'gtdlload'
require 'gthelper'
require 'libgtext/genome_stream'

module GT
  extend DL::Importable
  gtdlload "libgt"
  typealias "bool", "ibool"
  extern "GenomeStream* gff3_in_stream_new_sorted(const char *, bool)"

  class GFF3InStream
    include GenomeStream
    attr_reader :genome_stream
    def initialize(filename)
      if not File.readable?(filename)
        GT.gterror("file '#{filename}' not readable")
      end
      @genome_stream = GT.gff3_in_stream_new_sorted(filename, false)
      @genome_stream.free = GT::symbol("genome_stream_delete", "0P")
    end
  end
end
