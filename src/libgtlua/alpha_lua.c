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

#include "lauxlib.h"
#include "libgtcore/alpha.h"
#include "libgtext/luahelper.h"
#include "libgtlua/alpha_lua.h"

#define ALPHA_METATABLE  "GenomeTools.alpha"
#define check_alpha(L, POS) \
        (Alpha**) luaL_checkudata(L, POS, ALPHA_METATABLE)

static int alpha_lua_new_protein(lua_State *L)
{
  Alpha **alpha;
  assert(L);
  alpha = lua_newuserdata(L, sizeof (Alpha*));
  assert(alpha);
  *alpha = alpha_new_protein();
  assert(*alpha);
  luaL_getmetatable(L, ALPHA_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int alpha_lua_decode(lua_State *L)
{
  Alpha **alpha;
  unsigned int code;
  char character;
  alpha = check_alpha(L, 1);
  code = luaL_checkinteger(L, 2);
  /* XXX: too restrictive, does not consider wildcards */
  luaL_argcheck(L, code < alpha_size(*alpha), 2, "invalid code");
  character = alpha_decode(*alpha, code);
  lua_pushlstring(L, &character, 1);
  return 1;
}

static int alpha_lua_size(lua_State *L)
{
  Alpha **alpha;
  unsigned int size;
  alpha = check_alpha(L, 1);
  size = alpha_size(*alpha);
  lua_pushinteger(L, size);
  return 1;
}

static int alpha_lua_delete(lua_State *L)
{
  Alpha **alpha;
  alpha = check_alpha(L, 1);
  alpha_delete(*alpha);
  return 0;
}

static const struct luaL_Reg alpha_lib_f [] = {
  { "alpha_new_protein", alpha_lua_new_protein },
  { NULL, NULL }
};

static const struct luaL_Reg alpha_lib_m [] = {
  { "decode", alpha_lua_decode },
  { "size", alpha_lua_size },
  { NULL, NULL }
};

int luaopen_alpha(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, ALPHA_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, alpha_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, alpha_lib_m);
  luaL_register(L, "gt", alpha_lib_f);
  return 1;
}
