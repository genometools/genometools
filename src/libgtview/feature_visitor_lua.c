/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "lauxlib.h"
#include "gtlua.h"
#include "libgtext/genome_visitor_lua.h"
#include "libgtview/feature_index_lua.h"
#include "libgtview/feature_visitor.h"
#include "libgtview/feature_visitor_lua.h"

static int feature_visitor_lua_new(lua_State *L)
{
  GenomeVisitor **feature_visitor;
  FeatureIndex **feature_index;
  Env *env = get_env_from_registry(L);
  feature_visitor = lua_newuserdata(L, sizeof (GenomeVisitor**));
  assert(feature_visitor);
  feature_index = check_feature_index(L, 1);
  *feature_visitor = feature_visitor_new(*feature_index, env);
  luaL_getmetatable(L, GENOME_VISITOR_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static const struct luaL_Reg feature_visitor_lib_f [] = {
  { "feature_visitor_new", feature_visitor_lua_new },
  { NULL, NULL }
};

int luaopen_feature_visitor(lua_State *L)
{
  assert(L);
  luaL_register(L, "gt", feature_visitor_lib_f);
  return 1;
}
