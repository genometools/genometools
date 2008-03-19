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
  gtdlload "libgt"
  typealias "bool", "ibool"
  extern "FeatureIndex* feature_index_new()"
  extern "void feature_index_delete(FeatureIndex*)"
  extern "Array* feature_index_get_features_for_seqid(FeatureIndex*, const " +
                                                     "char*)"
  extern "const char* feature_index_get_first_seqid(const FeatureIndex*)"
  extern "StrArray* feature_index_get_seqids(const FeatureIndex*)"
  extern "void feature_index_get_rangeptr_for_seqid(FeatureIndex*, Range*, " +
                                                   "const char*)"
  extern "bool feature_index_has_seqid(const FeatureIndex*, const char*)"
  extern "void feature_index_delete(FeatureIndex*)"

  class FeatureIndex
    attr_reader :feature_index
    def initialize
      @feature_index = GT.feature_index_new()
      @feature_index.free = GT::symbol("feature_index_delete", "0P")
    end

    def get_features_for_seqid(seqid)
      rval = GT.feature_index_get_features_for_seqid(@feature_index, seqid)
      if rval then
        a = GT::Array.new(rval)
        result = []
        1.upto(a.size) do |i|
          result.push(GT::GenomeNode.new(GT.genome_node_rec_ref(a.get(i-1))))
        end
        result
      else
        nil
      end
    end

    def get_first_seqid
      GT.feature_index_get_first_seqid(@feature_index)
    end

    def get_seqids
      GT::StrArray.new(GT.feature_index_get_seqids(@feature_index)).to_a
    end

    def get_range_for_seqid(seqid)
      if not GT.feature_index_has_seqid(@feature_index, seqid)
        GT.gterror("feature_index does not contain seqid")
      end
      range = GT::Range.malloc
      GT.feature_index_get_rangeptr_for_seqid(@feature_index, range, seqid)
      range
    end
  end
end
