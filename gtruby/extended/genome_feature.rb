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

require 'dl/import'
require 'gthelper'
require 'core/str'
require 'extended/genome_node'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  typealias "bool", "ibool"
  extern "const GtGenomeNodeClass* gt_feature_node_class(void)"
  extern "void* gt_genome_node_cast(const GtGenomeNodeClass*, GtGenomeNode*)"
  extern "const char* gt_feature_node_get_type(GtFeatureNode*)"
  extern "int gt_feature_node_get_strand(GtFeatureNode*)"
  extern "int gt_feature_node_get_phase(GtFeatureNode*)"
  extern "float gt_feature_node_get_score(GtFeatureNode*)"
  extern "bool gt_feature_node_score_is_defined(const GtFeatureNode*)"
  extern "const char* gt_feature_node_get_attribute(GtGenomeNode*, " +
                                                     "const char*)"
  extern "void gt_feature_node_foreach_attribute(GtFeatureNode*, void*, " +
                                                  "void*)"

  #callback to populate attribute list
  def collect_attrib(tag, val, data)
    GT.gt_str_array_add_cstr(data, tag)
    GT.gt_str_array_add_cstr(data, val)
  end
  COLLECTFUNC = callback "void collect_attrib(const char*, const char*, void*)"

  class GenomeFeature < GenomeNode

    def initialize(gn, single=false)
      super(gn, single)
      attribs = GT::StrArray.new
      GT.gt_feature_node_foreach_attribute(@genome_node, COLLECTFUNC, attribs)
      attr_a = attribs.to_a
      @attribs = {}
      while not attr_a.empty? do
        @attribs[attr_a.shift] = attr_a.shift
      end
    end

    def get_type
      GT.gt_feature_node_get_type(@genome_node)
    end

    def get_strand
      GT.gt_feature_node_get_strand(@genome_node)
    end

    def get_phase
      GT.gt_feature_node_get_phase(@genome_node)
    end

    def get_score
      if GT.gt_feature_node_score_is_defined(@genome_node) then
        GT.gt_feature_node_get_score(@genome_node)
      else
        nil
      end
    end

    def get_attribute(attrib)
      @attribs[attrib]
    end

    def each_attribute
      @attribs.each_pair do |tag, val|
        yield tag, val
      end
    end

    def to_ptr
      @genome_node
    end
  end
end
