/*
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
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
#include "extended/feature_index_api.h"

/* exports the FeatureIndex class to Lua:

   -- Returns a new FeatureIndex object storing the index in memory.
   function feature_index_memory_new()

   -- Add all features from all sequence regions contained in <gff3file> to
   -- <feature_index>.
   function feature_index:add_gff3file(gff3file)

   -- Add <region_node> to <feature_index>.
   function feature_index:add_region_node(region_node)

   -- Add <feature_node> to <feature_index>, implicitly creating sequence
   -- region if not present before.
   function feature_index:add_feature_node(feature_node)

   -- Returns the feature nodes for sequence ID <seqid> in an array.
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
int gt_lua_open_feature_index(lua_State*);

/* Push a <GtFeatureIndex*> to Lua, takes ownership! */
void gt_lua_feature_index_push(lua_State *L, GtFeatureIndex *fi);

#define FEATURE_INDEX_METATABLE  "GenomeTools.feature_index"
#define check_feature_index(L, POS) \
          (GtFeatureIndex**) luaL_checkudata(L, POS, FEATURE_INDEX_METATABLE)

#endif
