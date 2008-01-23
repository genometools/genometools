/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef FEATURE_INDEX_LUA_H
#define FEATURE_INDEX_LUA_H

#include "lua.h"

/* exports the FeatureIndex class to Lua:

   -- Returns a new <feature_index> object.
   function feature_index_new()

   -- Add <sequence_region> to <feature_index>.
   function feature_index:add_sequence_region(sequence_region)

   -- Add <genome_feature> to <feature_index>.
   function feature_index:add_genome_feature(genome_feature)

   -- Returns the genome features for sequence ID <seqid> in an array.
   function feature_index:get_features_for_seqid(seqid)

   -- Returns the genome features for sequence ID <seqid> within <range> in an
   -- array.
   function feature_index:get_features_for_range(seqid, range)

   -- Returns the first sequence ID stored in <feature_index>.
   function feature_index:get_first_seqid()

   -- Returns an array containins all sequence IDs stored in <feature_index>.
   function feature_index:get_seqids()

   -- Returns the range covered by features of sequence ID <seqid> in
   -- <feature_index>.
   function feature_index:get_range_for_seqid(seqid)
*/
int luaopen_feature_index(lua_State*);

#define FEATURE_INDEX_METATABLE  "GenomeTools.feature_index"
#define check_feature_index(L, POS) \
          (FeatureIndex**) luaL_checkudata(L, POS, FEATURE_INDEX_METATABLE)

#endif
