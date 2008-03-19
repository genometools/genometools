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
require 'libgtcore/range'

module GT
  extend DL::Importable
  gtdlload "libgt"
  extern "Diagram* diagram_new(FeatureIndex*, const char*, const Range*, " +
                              "Config*)"
  extern "void diagram_delete(Diagram*)"

  class Diagram
    attr_reader :diagram
    def initialize(feature_index, seqid, range, config)
      if range.start > range.end
        GT.gterror("range.start > range.end")
      end
      @diagram = GT.diagram_new(feature_index.feature_index, seqid, range,
                                config.config)
      @diagram.free = GT::symbol("diagram_delete", "0P")
    end
  end
end
