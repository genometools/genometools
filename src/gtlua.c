/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "gtlua.h"
#include "lauxlib.h"
#include "libgtcore/gtcorelua.h"
#include "libgtext/gtextlua.h"

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

int luaopen_gt(lua_State *L)
{
  assert(L);
  luaopen_gtcore(L); /* open core library */
  luaopen_gtext(L);  /* open extended library */
  return 1;
}

void run_interactive_lua_interpreter(lua_State *L)
{
  char buf[BUFSIZ];
  int error;
  assert(L);
  while (fgets(buf, sizeof buf, stdin)) {
    error = luaL_loadbuffer(L, buf, strlen(buf), "line") ||
            lua_pcall(L, 0, 0, 0);
    if (error) {
      fprintf(stderr, "%s", lua_tostring(L, -1));
      lua_pop(L, 1); /* pop error message */
    }
  }
}
