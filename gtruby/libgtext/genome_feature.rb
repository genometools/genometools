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
require 'libgtcore/str'
require 'libgtext/genome_node'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  typealias "bool", "ibool"
  extern "const GT_GenomeNodeClass* gt_genome_feature_class(void)"
  extern "void* gt_genome_node_cast(const GT_GenomeNodeClass*, GT_GenomeNode*)"
  extern "GT_FeatureType* gt_genome_feature_get_type(GT_GenomeFeature*)"
  extern "int gt_genome_feature_get_strand(GT_GenomeFeature*)"
  extern "int gt_genome_feature_get_phase(GT_GenomeFeature*)"
  extern "float gt_genome_feature_get_score(GT_GenomeFeature*)"
  extern "bool gt_genome_feature_score_is_defined(const GT_GenomeFeature*)"
  extern "const char* gt_feature_type_get_cstr(const GT_FeatureType*)"
  extern "const char* gt_genome_feature_get_attribute(GT_GenomeNode*,const char*)"
  extern "void gt_genome_feature_foreach_attribute(GT_GenomeFeature*, void*, void*)"

  #callback to populate attribute list
  def collect_attrib(tag, val, data)
    GT.gt_strarray_add_cstr(data, tag)
    GT.gt_strarray_add_cstr(data, val)
  end
  COLLECTFUNC = callback "void collect_attrib(const char*, const char*, void*)"

  class GenomeFeature < GenomeNode

    def initialize(gn, single=false)
      super(gn, single)
      attribs = GT::StrArray.new
      GT.gt_genome_feature_foreach_attribute(@genome_node, COLLECTFUNC, attribs)
      attr_a = attribs.to_a
      @attribs = {}
      while not attr_a.empty? do
        @attribs[attr_a.shift] = attr_a.shift
      end
    end

    def get_type
      GT.gt_feature_type_get_cstr(GT.gt_genome_feature_get_type(@genome_node))
    end

    def get_strand
      GT.gt_genome_feature_get_strand(@genome_node)
    end

    def get_phase
      GT.gt_genome_feature_get_phase(@genome_node)
    end

    def get_score
      if GT.gt_genome_feature_score_is_defined(@genome_node) then
        GT.gt_genome_feature_get_score(@genome_node)
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
