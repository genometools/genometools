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
  gtdlload "libgt"
  extern "const GenomeNodeClass* genome_feature_class(void)"
  extern "void* genome_node_cast(const GenomeNodeClass*, GenomeNode*)"
  extern "GenomeFeatureType* genome_feature_get_type(GenomeFeature*)"
  extern "int genome_feature_get_strand(GenomeFeature*)"
  extern "int genome_feature_get_phase(GenomeFeature*)"
  extern "double genome_feature_get_score(GenomeFeature*)"
  extern "const char* genome_feature_type_get_cstr(const GenomeFeatureType*)"
  extern "const char* genome_feature_get_attribute(GenomeNode*,const char*)"
  extern "void genome_feature_foreach_attribute(GenomeFeature*, void*, void*)"

  #callback to populate attribute list
  def collect_attrib(tag, val, data)
    GT.strarray_add_cstr(data, tag)
    GT.strarray_add_cstr(data, val)
  end
  COLLECTFUNC = callback "void collect_attrib(const char*, const char*, void*)"

  class GenomeFeature < GenomeNode

    def initialize(gn, single=false)
      super(gn, single)
      attribs = GT::StrArray.new
      GT.genome_feature_foreach_attribute(@genome_node, COLLECTFUNC, attribs)
      attr_a = attribs.to_a
      @attribs = {}
      while not attr_a.empty? do
        @attribs[attr_a.shift] = attr_a.shift
      end
    end

    def get_type
      GT.genome_feature_type_get_cstr(GT.genome_feature_get_type(@genome_node))
    end

    def get_strand
      GT.genome_feature_get_strand(@genome_node)
    end

    def get_phase
      GT.genome_feature_get_phase(@genome_node)
    end

    def get_score
      GT.genome_feature_get_score(@genome_node)
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
