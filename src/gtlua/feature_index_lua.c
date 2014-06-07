/*
  Copyright (c) 2007 Gordon Gremme <gordon@gremme.org>
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
#include "core/unused_api.h"
#include "core/ma.h"
#include "extended/feature_index_api.h"
#include "extended/feature_index_memory_api.h"
#include "extended/feature_node.h"
#include "extended/luahelper.h"
#include "gtlua/feature_index_lua.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/range_lua.h"

static int feature_index_memory_lua_new(lua_State *L)
{
  GtFeatureIndex **feature_index;
  feature_index = lua_newuserdata(L, sizeof (GtFeatureIndex*));
  gt_assert(feature_index);
  *feature_index = gt_feature_index_memory_new();
  luaL_getmetatable(L, FEATURE_INDEX_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int feature_index_lua_add_region_node(lua_State *L)
{
  GtFeatureIndex **fi;
  GtGenomeNode **gn;
  GtRegionNode *rn;
  GtError *err;
  gt_assert(L);
  fi = check_feature_index(L, 1);
  gn = check_genome_node(L, 2);
  rn = gt_region_node_try_cast(*gn);
  luaL_argcheck(L, rn, 2, "not a region node");
  err = gt_error_new();
  if (gt_feature_index_add_region_node(*fi, rn, err))
    return gt_lua_error(L, err);
  gt_error_delete(err);
  return 0;
}

static int feature_index_lua_add_gff3file(lua_State *L)
{
  GtFeatureIndex **fi;
  const char *filename;
  GtError *err;
  gt_assert(L);
  fi = check_feature_index(L, 1);
  filename = luaL_checkstring(L, 2);
  err = gt_error_new();
  if (gt_feature_index_add_gff3file(*fi, filename, err))
    return gt_lua_error(L, err);
  gt_error_delete(err);
  return 0;
}

static int feature_index_lua_add_feature_node(lua_State *L)
{
  GtFeatureIndex **fi;
  GtGenomeNode **gn;
  GtFeatureNode *fn;
  GtStr *seqid;
  GtError *err;
  bool has_seqid;
  gt_assert(L);
  fi = check_feature_index(L, 1);
  gn = check_genome_node(L, 2);
  fn = gt_feature_node_cast(*gn);
  luaL_argcheck(L, fn, 2, "not a feature node");
  seqid = gt_genome_node_get_seqid(*gn);
  luaL_argcheck(L, seqid, 2, "feature does not have a sequence id");
  err = gt_error_new();
  if (gt_feature_index_has_seqid(*fi, &has_seqid, gt_str_get(seqid), err))
    return gt_lua_error(L, err);
  gt_error_delete(err);
  luaL_argcheck(L, has_seqid, 2,
                "feature index does not contain corresponding sequence region");
  err = gt_error_new();
  if (gt_feature_index_add_feature_node(*fi, fn, err))
    return gt_lua_error(L, err);
  gt_error_delete(err);
  return 0;
}

static void push_features_as_table(lua_State *L, GtArray *features)
{
  GtUword i;
  if (features && gt_array_size(features)) {
    /* push table containing feature references onto the stack */
    lua_newtable(L);
    for (i = 0; i < gt_array_size(features); i++) {
      lua_pushinteger(L, i+1); /* in Lua we index from 1 on */
      gt_lua_genome_node_push(L, gt_genome_node_ref(*(GtGenomeNode**)
                                                  gt_array_get(features, i)));
      lua_rawset(L, -3);
    }
  }
  else
    lua_pushnil(L);
}

static int feature_index_lua_get_features_for_seqid(lua_State *L)
{
  GtFeatureIndex **feature_index;
  const char *seqid;
  GtArray *features;
  GtError *err;
  feature_index = check_feature_index(L, 1);
  seqid = luaL_checkstring(L, 2);
  err = gt_error_new();
  features = gt_feature_index_get_features_for_seqid(*feature_index, seqid,
                                                     err);
  if (!features)
    return gt_lua_error(L, err);
  gt_error_delete(err);
  push_features_as_table(L, features);
  gt_array_delete(features);
  return 1;
}

static int feature_index_lua_get_features_for_range(lua_State *L)
{
  GtFeatureIndex **feature_index;
  const char *seqid;
  GtRange *range;
  GtError *err;
  bool has_seqid;
  GtArray *features;
  GT_UNUSED int had_err;

  feature_index = check_feature_index(L, 1);
  seqid = luaL_checkstring(L, 2);
  err = gt_error_new();
  if (gt_feature_index_has_seqid(*feature_index, &has_seqid, seqid, err))
    return gt_lua_error(L, err);
  gt_error_delete(err);
  luaL_argcheck(L, has_seqid, 2,
                "feature_index does not contain seqid");
  range = check_range(L, 3);
  features = gt_array_new(sizeof (GtGenomeNode*));
  err = gt_error_new();
  had_err = gt_feature_index_get_features_for_range(*feature_index, features,
                                                    seqid, range, err);
  if (had_err)
    return gt_lua_error(L, err);
  gt_error_delete(err);
  push_features_as_table(L, features);
  gt_array_delete(features);
  return 1;
}

static int feature_index_lua_get_first_seqid(lua_State *L)
{
  GtFeatureIndex **feature_index;
  char *seqid;
  GtError *err;
  feature_index = check_feature_index(L, 1);
  err = gt_error_new();
  seqid = gt_feature_index_get_first_seqid(*feature_index, err);
  if (gt_error_is_set(err))
    return gt_lua_error(L, err);
  gt_error_delete(err);
  if (seqid) {
    lua_pushstring(L, seqid);
    gt_free(seqid);
  } else
    lua_pushnil(L);
  return 1;
}

static int feature_index_lua_get_seqids(lua_State *L)
{
  GtFeatureIndex **feature_index;
  GtStrArray *seqids;
  GtError *err;
  feature_index = check_feature_index(L, 1);
  err = gt_error_new();
  seqids = gt_feature_index_get_seqids(*feature_index, err);
  if (!seqids)
    return gt_lua_error(L, err);
  gt_error_delete(err);
  /* push table containing sequence ids onto the stack */
  gt_lua_push_strarray_as_table(L, seqids);
  gt_str_array_delete(seqids);
  return 1;
}

static int feature_index_lua_get_range_for_seqid(lua_State *L)
{
  GtFeatureIndex **feature_index;
  const char *seqid;
  GtError *err;
  GtRange range;
  bool has_seqid;
  feature_index = check_feature_index(L, 1);
  seqid = luaL_checkstring(L, 2);
  err = gt_error_new();
  if (gt_feature_index_has_seqid(*feature_index, &has_seqid, seqid, err))
    return gt_lua_error(L, err);
  gt_error_delete(err);
  luaL_argcheck(L, has_seqid, 2,
                "feature_index does not contain seqid");
  err = gt_error_new();
  if (gt_feature_index_get_range_for_seqid(*feature_index, &range, seqid, err))
    return gt_lua_error(L, err);
  gt_error_delete(err);
  return gt_lua_range_push(L, range);
}

void gt_lua_feature_index_push(lua_State *L, GtFeatureIndex *fi)
{
  GtFeatureIndex **fi_lua;
  gt_assert(L && fi);
  fi_lua = lua_newuserdata(L, sizeof (GtFeatureIndex**));
  *fi_lua = fi;
  luaL_getmetatable(L, FEATURE_INDEX_METATABLE);
  lua_setmetatable(L, -2);
}

static int feature_index_lua_delete(lua_State *L)
{
  GtFeatureIndex **feature_index = check_feature_index(L, 1);
  gt_feature_index_delete(*feature_index);
  return 0;
}

static const struct luaL_Reg feature_index_lib_f [] = {
  { "feature_index_memory_new", feature_index_memory_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg feature_index_lib_m [] = {
  { "add_region_node", feature_index_lua_add_region_node },
  { "add_feature_node", feature_index_lua_add_feature_node },
  { "add_gff3file", feature_index_lua_add_gff3file },
  { "get_features_for_seqid", feature_index_lua_get_features_for_seqid },
  { "get_features_for_range", feature_index_lua_get_features_for_range },
  { "get_first_seqid", feature_index_lua_get_first_seqid },
  { "get_seqids", feature_index_lua_get_seqids },
  { "get_range_for_seqid", feature_index_lua_get_range_for_seqid },
  { NULL, NULL }
};

int gt_lua_open_feature_index(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
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
  gt_lua_export_metatable(L, FEATURE_INDEX_METATABLE);
  luaL_register(L, "gt", feature_index_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}
