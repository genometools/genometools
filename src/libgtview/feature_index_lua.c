/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "lauxlib.h"
#include "gtlua.h"
#include "libgtext/genome_node_lua.h"
#include "libgtview/feature_index.h"
#include "libgtview/feature_index_lua.h"

static int feature_index_lua_new(lua_State *L)
{
  FeatureIndex **feature_index;
  Env *env = get_env_from_registry(L);
  feature_index = lua_newuserdata(L, sizeof (FeatureIndex**));
  assert(feature_index);
  *feature_index = feature_index_new(env);
  luaL_getmetatable(L, FEATURE_INDEX_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static void push_features_as_table(lua_State *L, Array *features, Env *env)
{
  unsigned long i;
  if (features && array_size(features)) {
    /* push table containing feature references onto the stack */
    lua_newtable(L);
    for (i = 0; i < array_size(features); i++) {
      lua_pushinteger(L, i+1); /* in Lua we index from 1 on */
      genome_node_lua_push(L, genome_node_rec_ref(*(GenomeNode**)
                                                  array_get(features, i), env));
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
  Env *env = get_env_from_registry(L);
  feature_index = check_feature_index(L, 1);
  seqid = luaL_checkstring(L, 2);
  features = feature_index_get_features_for_seqid(*feature_index, seqid);
  push_features_as_table(L, features, env);
  return 1;
}

static int feature_index_lua_get_features_for_range(lua_State *L)
{
  FeatureIndex **feature_index;
  const char *seqid;
  Range range;
  Array *features;
  int had_err;
  Env *env = get_env_from_registry(L);
  feature_index = check_feature_index(L, 1);
  seqid = luaL_checkstring(L, 2);
  luaL_argcheck(L, feature_index_has_seqid(*feature_index, seqid, env), 2,
                "feature_index does not contain seqid");
  range.start = luaL_checklong(L, 3);
  range.end   = luaL_checklong(L, 4);
  luaL_argcheck(L, range.start <= range.end, 3, "must be <= endpos");
  features = array_new(sizeof (GenomeNode*), env);
  had_err = feature_index_get_features_for_range(*feature_index, features,
                                                 seqid, range, env);
  assert(!had_err); /* it was checked before that the feature_index contains the
                       given sequence id*/
  push_features_as_table(L, features, env);
  array_delete(features, env);
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

static int feature_index_lua_get_range_for_seqid(lua_State *L)
{
  FeatureIndex **feature_index;
  const char *seqid;
  Range range;
  Env *env = get_env_from_registry(L);
  feature_index = check_feature_index(L, 1);
  seqid = luaL_checkstring(L, 2);
  luaL_argcheck(L, feature_index_has_seqid(*feature_index, seqid, env), 2,
                "feature_index does not contain seqid");
  range = feature_index_get_range_for_seqid(*feature_index, seqid);
  lua_pushinteger(L, range.start);
  lua_pushinteger(L, range.end);
  return 2;
}

static int feature_index_lua_delete(lua_State *L)
{
  FeatureIndex **feature_index = check_feature_index(L, 1);
  Env *env;
  env = get_env_from_registry(L);
  feature_index_delete(*feature_index, env);
  return 0;
}

static const struct luaL_Reg feature_index_lib_f [] = {
  { "feature_index_new", feature_index_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg feature_index_lib_m [] = {
  { "get_features_for_seqid", feature_index_lua_get_features_for_seqid },
  { "get_features_for_range", feature_index_lua_get_features_for_range },
  { "get_first_seqid", feature_index_lua_get_first_seqid },
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
  luaL_register(L, "gt", feature_index_lib_f);
  return 1;
}
