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

module GT
  extend DL::Importable
  gtdlload "libgt"
  extern "int genome_node_accept(GenomeNode*, GenomeVisitor*, Error*)"
  extern "GenomeNode* genome_node_rec_ref(GenomeNode*)"
  extern "void genome_node_rec_delete(GenomeNode*)"

  class GenomeNode
    attr_reader :genome_node
    def initialize(node_ptr)
      @genome_node = node_ptr
      @genome_node.free = GT::symbol("genome_node_rec_delete", "0P")
    end

    def accept(visitor)
      err = GT::Error.new()
      rval = GT.genome_node_accept(self.genome_node, visitor.genome_visitor,
                                   err.to_ptr)
      if rval != 0
        GT.gterror(err)
      end
    end
  end
end
