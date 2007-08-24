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
#include "libgtcore/gtcore_lua.h"
#include "libgtext/gtext_lua.h"

#ifdef LIBGTVIEW
#include "libgtview/gtview_lua.h"
#endif

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
#ifdef LIBGTVIEW
  luaopen_gtview(L); /* open view library */
#endif
  return 1;
}

void set_arg_in_lua_interpreter(lua_State *L, const char *argv_0,
                                const char **argv)
{
  lua_Integer n = 0;
  assert(L && argv_0);
  /* create table */
  lua_newtable(L);
  /* set arg[0] */
  lua_pushinteger(L, 0);
  lua_pushstring(L, argv_0);
  lua_rawset(L, -3);
  /* set other arguments */
  while (argv[n]) {
    lua_pushinteger(L, n+1);
    lua_pushstring(L, argv[n]);
    lua_rawset(L, -3);
    n++;
  }
  /* register table globally */
  lua_setglobal(L, "arg");
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

int luagt_error(lua_State *L, Env *env)
{
  assert(L && env);
  assert(env_error_is_set(env));
  luaL_where(L, 1);
  lua_pushstring(L, env_error_get(env));
  env_error_unset(env);
  lua_concat(L, 2);
  return lua_error(L);
}
