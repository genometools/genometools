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

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "lauxlib.h"

#ifdef LIBGTVIEW
#include "libgtview/config.h"
#endif

/* key used to store the Env object in the Lua registry */
#define ENV_KEY env_new

#ifdef LIBGTVIEW
/* key used to store the Config object in the Lua registry */
#define CONFIG_KEY config_new

void put_config_in_registry(lua_State *L, Config *config)
{
  assert(L && config);
  lua_pushlightuserdata(L, CONFIG_KEY);
  lua_pushlightuserdata(L, config);
  lua_rawset(L, LUA_REGISTRYINDEX);
}

Config* get_config_from_registry(lua_State *L)
{
  Config *config;
  assert(L);
  lua_pushlightuserdata(L, CONFIG_KEY);
  lua_rawget(L, LUA_REGISTRYINDEX);
  assert(lua_islightuserdata(L, -1));
  config = lua_touserdata(L, -1);
  lua_pop(L, 1);
  return config;
}
#endif

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
  lua_pop(L, 1);
  return env;
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
      fprintf(stderr, "%s\n", lua_tostring(L, -1));
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
