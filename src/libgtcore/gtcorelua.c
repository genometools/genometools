/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "libgtcore/bittablua.h"
#include "libgtcore/gtcorelua.h"

/* key used to store the Env object in the Lua registry */
#define ENV_KEY env_new

void put_env_in_registry(lua_State *L, Env *env)
{
  assert(L && env);
  lua_pushlightuserdata(L, ENV_KEY); /* push the key */
  lua_pushlightuserdata(L, env); /* push the value */
  lua_rawset(L, LUA_REGISTRYINDEX); /* store env in registry */
}

Env* get_env_from_registry(lua_State *L)
{
  Env *env;
  assert(L);
  lua_pushlightuserdata(L, ENV_KEY);
  lua_rawget(L, LUA_REGISTRYINDEX);
  assert(lua_islightuserdata(L, -1));
  env = lua_touserdata(L, -1);
  assert(env);
  lua_pop(L, 1);
  return env;
}

int luaopen_gtcore(lua_State *L)
{
  assert(L);
  luaopen_bittab(L);
  return 1;
}
