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

#define FEATURE_INDEX_METATABLE  "GenomeTools.feature_index"
#define checkfeature_index(L) \
        (FeatureIndex**) luaL_checkudata(L, 1, FEATURE_INDEX_METATABLE);

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

static int feature_index_lua_get_features_for_seqid(lua_State *L)
{
  FeatureIndex **feature_index;
  const char *seqid;
  Array *features;
  unsigned long i;
  Env *env = get_env_from_registry(L);
  feature_index = checkfeature_index(L);
  seqid = luaL_checkstring(L, 2);
  features = feature_index_get_features_for_seqid(*feature_index, seqid);
  if (features) {
    /* push table containing feature references onto the stack */
    assert(array_size(features));
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
  return 1;
}

static int feature_index_lua_delete(lua_State *L)
{
  FeatureIndex **feature_index = checkfeature_index(L);
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
