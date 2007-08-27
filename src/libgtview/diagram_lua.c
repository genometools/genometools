/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "lauxlib.h"
#include "gtlua.h"
#include "libgtext/genome_node_lua.h"
#include "libgtview/diagram.h"
#include "libgtview/diagram_lua.h"

static int diagram_lua_new(lua_State *L)
{
  Diagram **diagram;
  lua_Integer i = 1;
  Array *features;
  GenomeNode **gn;
  Range range;
  const char *msg;
  bool error;
  Config *config;
  Env *env = get_env_from_registry(L);
  diagram = lua_newuserdata(L, sizeof (Diagram**));
  assert(diagram);
  /* make sure we got a table as first argument */
  luaL_checktype(L, 1, LUA_TTABLE);
  /* traverse table and save the features */
  features = array_new(sizeof (GenomeNode*), env);
  lua_pushinteger(L, i);
  lua_gettable(L, 1);
  while (!lua_isnil(L, -1)) {
    error = false;
    gn = lua_touserdata(L, -1);
    if (gn && lua_getmetatable(L, -1)) {
      lua_getfield(L, LUA_REGISTRYINDEX, GENOME_NODE_METATABLE);
      if (lua_rawequal(L, -1, -2)) {
        lua_pop(L, 2); /* remove both metatables */
        array_add(features, *gn, env);
      }
      else
        error = true;
    }
    else
      error = true;
    if (error) {
      /* we do not have a genome node in the table */
      msg = lua_pushfstring(L, "expected %s as type of feature number %d",
                            GENOME_NODE_METATABLE, i);
      array_delete(features, env);
      return lua_error(L);
    }
    i++;
    lua_pushinteger(L, i);
    lua_gettable(L, 1);
  }
  /* get range */
  range.start = luaL_checklong(L, 2);
  range.end   = luaL_checklong(L, 3);
  luaL_argcheck(L, range.start <= range.end, 2, "must be <= endpos");
  /* create diagram */
  config = get_config_from_registry(L);
  *diagram = diagram_new(features, range, config, env);
  array_delete(features, env);
  luaL_getmetatable(L, DIAGRAM_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int diagram_lua_delete(lua_State *L)
{
  Diagram **diagram;
  Env *env;
  diagram = check_diagram(L, 1);
  env = get_env_from_registry(L);
  diagram_delete(*diagram, env);
  return 0;
}

static const struct luaL_Reg diagram_lib_f [] = {
  { "diagram_new", diagram_lua_new },
  { NULL, NULL }
};

int luaopen_diagram(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, DIAGRAM_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, diagram_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, "gt", diagram_lib_f);
  return 1;
}
