/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "lauxlib.h"
#include "libgtext/luahelper.h"
#include "libgtlua/feature_index_lua.h"
#include "libgtlua/genome_node_lua.h"
#include "libgtlua/range_lua.h"
#include "libgtview/feature_index.h"

static int feature_index_lua_new(lua_State *L)
{
  FeatureIndex **feature_index;
  feature_index = lua_newuserdata(L, sizeof (FeatureIndex*));
  assert(feature_index);
  *feature_index = feature_index_new();
  luaL_getmetatable(L, FEATURE_INDEX_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int feature_index_lua_add_sequence_region(lua_State *L)
{
  FeatureIndex **fi;
  GenomeNode **gn;
  SequenceRegion *sr;
  assert(L);
  fi = check_feature_index(L, 1);
  gn = check_genome_node(L, 2);
  sr = genome_node_cast(sequence_region_class(), *gn);
  luaL_argcheck(L, sr, 2, "not a sequence region");
  feature_index_add_sequence_region(*fi, sr);
  return 0;
}

static int feature_index_lua_add_genome_feature(lua_State *L)
{
  FeatureIndex **fi;
  GenomeNode **gn;
  GenomeFeature *gf;
  Str *seqid;
  assert(L);
  fi = check_feature_index(L, 1);
  gn = check_genome_node(L, 2);
  gf = genome_node_cast(genome_feature_class(), *gn);
  luaL_argcheck(L, gf, 2, "not a genome feature");
  seqid = genome_node_get_seqid(*gn);
  luaL_argcheck(L, seqid, 2, "genome_feature does not have a sequence id");
  luaL_argcheck(L, feature_index_has_seqid(*fi, str_get(seqid)), 2,
                "feature index does not contain corresponding sequence region");
  feature_index_add_genome_feature(*fi, gf);
  return 0;
}

static void push_features_as_table(lua_State *L, Array *features)
{
  unsigned long i;
  if (features && array_size(features)) {
    /* push table containing feature references onto the stack */
    lua_newtable(L);
    for (i = 0; i < array_size(features); i++) {
      lua_pushinteger(L, i+1); /* in Lua we index from 1 on */
      genome_node_lua_push(L, genome_node_rec_ref(*(GenomeNode**)
                                                  array_get(features, i)));
      lua_rawset(L, -3);
    }
  }
  else
    lua_pushnil(L);
}

static int feature_index_lua_get_features_for_seqid(lua_State *L)
{
  FeatureIndex **feature_index;
  const char *seqid;
  Array *features;
  feature_index = check_feature_index(L, 1);
  seqid = luaL_checkstring(L, 2);
  features = feature_index_get_features_for_seqid(*feature_index, seqid);
  push_features_as_table(L, features);
  array_delete(features);
  return 1;
}

static int feature_index_lua_get_features_for_range(lua_State *L)
{
  FeatureIndex **feature_index;
  const char *seqid;
  Range *range;
  Array *features;
  int had_err;
  feature_index = check_feature_index(L, 1);
  seqid = luaL_checkstring(L, 2);
  luaL_argcheck(L, feature_index_has_seqid(*feature_index, seqid), 2,
                "feature_index does not contain seqid");
  range = check_range(L, 3);
  features = array_new(sizeof (GenomeNode*));
  had_err = feature_index_get_features_for_range(*feature_index, features,
                                                 seqid, *range, NULL);
  assert(!had_err); /* it was checked before that the feature_index contains the
                       given sequence id*/
  push_features_as_table(L, features);
  array_delete(features);
  return 1;
}

static int feature_index_lua_get_first_seqid(lua_State *L)
{
  FeatureIndex **feature_index;
  const char *seqid;
  feature_index = check_feature_index(L, 1);
  seqid = feature_index_get_first_seqid(*feature_index);
  if (seqid)
    lua_pushstring(L, seqid);
  else
    lua_pushnil(L);
  return 1;
}

static int feature_index_lua_get_seqids(lua_State *L)
{
  FeatureIndex **feature_index;
  StrArray *seqids;
  feature_index = check_feature_index(L, 1);
  seqids = feature_index_get_seqids(*feature_index);
  assert(seqids);
  /* push table containing sequence ids onto the stack */
  lua_push_strarray_as_table(L, seqids);
  strarray_delete(seqids);
  return 1;
}

static int feature_index_lua_get_range_for_seqid(lua_State *L)
{
  FeatureIndex **feature_index;
  const char *seqid;
  feature_index = check_feature_index(L, 1);
  seqid = luaL_checkstring(L, 2);
  luaL_argcheck(L, feature_index_has_seqid(*feature_index, seqid), 2,
                "feature_index does not contain seqid");
  return range_lua_push(L, feature_index_get_range_for_seqid(*feature_index,
                                                             seqid));
}

static int feature_index_lua_delete(lua_State *L)
{
  FeatureIndex **feature_index = check_feature_index(L, 1);
  feature_index_delete(*feature_index);
  return 0;
}

static const struct luaL_Reg feature_index_lib_f [] = {
  { "feature_index_new", feature_index_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg feature_index_lib_m [] = {
  { "add_sequence_region", feature_index_lua_add_sequence_region },
  { "add_genome_feature", feature_index_lua_add_genome_feature },
  { "get_features_for_seqid", feature_index_lua_get_features_for_seqid },
  { "get_features_for_range", feature_index_lua_get_features_for_range },
  { "get_first_seqid", feature_index_lua_get_first_seqid },
  { "get_seqids", feature_index_lua_get_seqids },
  { "get_range_for_seqid", feature_index_lua_get_range_for_seqid },
  { NULL, NULL }
};

int luaopen_feature_index(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, FEATURE_INDEX_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, feature_index_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, feature_index_lib_m);
  lua_export_metatable(L, FEATURE_INDEX_METATABLE);
  luaL_register(L, "gt", feature_index_lib_f);
  return 1;
}
