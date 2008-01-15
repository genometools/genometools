/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtview/config.h"

/* key used to store the Config object in the Lua registry */
#define CONFIG_KEY config_new

void lua_put_config_in_registry(lua_State *L, Config *config)
{
  assert(L && config);
  lua_pushlightuserdata(L, CONFIG_KEY);
  lua_pushlightuserdata(L, config);
  lua_rawset(L, LUA_REGISTRYINDEX);
}

Config* lua_get_config_from_registry(lua_State *L)
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
