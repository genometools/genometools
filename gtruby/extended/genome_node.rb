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
  gtdlload "libgenometools"
  extern "int gt_genome_node_accept(GT_GenomeNode*, GenomeVisitor*, Error*)"
  extern "GT_GenomeNode* gt_genome_node_rec_ref(GT_GenomeNode*)"
  extern "GT_GenomeNode* gt_genome_node_ref(GT_GenomeNode*)"
  extern "unsigned long gt_genome_node_get_start(GT_GenomeNode*)"
  extern "unsigned long gt_genome_node_get_end(GT_GenomeNode*)"
  extern "const char* gt_genome_node_get_filename(GT_GenomeNode*)"
  extern "void gt_genome_node_rec_delete(GT_GenomeNode*)"
  extern "void gt_genome_node_delete(GT_GenomeNode*)"

  class GenomeNode
    attr_reader :genome_node
    def initialize(node_ptr, single=false)
      # use 'single' if not referencing root nodes
      if single then
        @genome_node = GT.gt_genome_node_ref(node_ptr)
        @genome_node.free = GT::symbol("gt_genome_node_delete", "0P")
      else
        @genome_node = node_ptr
        @genome_node.free = GT::symbol("gt_genome_node_rec_delete", "0P")
      end
    end

    def get_range
      (GT::gt_genome_node_get_start(@genome_node)..\
        GT::gt_genome_node_get_end(@genome_node))
    end

    def get_filename
      GT.gt_genome_node_get_filename(@genome_node)
    end

    def to_ptr
      @genome_node
    end

    def accept(visitor)
      err = GT::Error.new()
      rval = GT.gt_genome_node_accept(@genome_node, visitor.genome_visitor,
                                      err.to_ptr)
      if rval != 0
        GT.gterror(err)
      end
    end
  end
end
