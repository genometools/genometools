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

require 'dl/import'
require 'libgtcore/array'
require 'libgtcore/range'
require 'libgtcore/strarray'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  typealias "bool", "ibool"
  extern "GtFeatureIndex* gt_feature_index_memory_new()"
  extern "void gt_feature_index_delete(GtFeatureIndex*)"
  extern "int gt_feature_index_add_gff3file(GtFeatureIndex*, " +
                                           "const char*, " +
                                           "GtError*)"
  extern "GtArray* gt_feature_index_get_features_for_seqid(Gt_FeatureIndex*, " +
                                                          "const char*)"
  extern "const char* gt_feature_index_get_first_seqid(const GtFeatureIndex*)"
  extern "GtStrArray* gt_feature_index_get_seqids(const GtFeatureIndex*)"
  extern "void gt_feature_index_get_range_for_seqid(GtFeatureIndex*, " +
                                                   "GtRange*, const char*)"
  extern "bool gt_feature_index_has_seqid(const GtFeatureIndex*, const char*)"
  extern "void gt_feature_index_delete(GtFeatureIndex*)"

  class FeatureIndex
    attr_reader :feature_index
    def initialize
      @feature_index = GT.gt_feature_index_memory_new()
      @feature_index.free = GT::symbol("gt_feature_index_delete", "0P")
    end

    def get_features_for_seqid(seqid)
      rval = GT.gt_feature_index_get_features_for_seqid(@feature_index, seqid)
      if rval then
        a = GT::Array.new(rval)
        result = []
        1.upto(a.size) do |i|
          result.push(GT::GenomeNode.new(GT.gt_genome_node_rec_ref(a.get(i-1))))
        end
        result
      else
        nil
      end
    end

    def add_gff3file(filename)
      err = GT::Error.new()
      rval = GT.gt_feature_index_add_gff3file(@feature_index, filename, \
                                              err.to_ptr)
      if rval != 0 then GT.gterror(err) end
    end

    def get_first_seqid
      GT.gt_feature_index_get_first_seqid(@feature_index)
    end

    def get_seqids
      GT::StrArray.new(GT.gt_feature_index_get_seqids(@feature_index)).to_a
    end

    def get_range_for_seqid(seqid)
      if not GT.gt_feature_index_has_seqid(@feature_index, seqid)
        GT.gterror("feature_index does not contain seqid")
      end
      range = GT::Range.malloc
      GT.gt_feature_index_get_range_for_seqid(@feature_index, range, seqid)
      range
    end
  end
end
