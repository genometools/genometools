 /*
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <string.h>
#include "lauxlib.h"
#include "core/alphabet.h"
#include "extended/luahelper.h"
#include "gtlua/alphabet_lua.h"

static int alphabet_lua_new_protein(lua_State *L)
{
  GtAlphabet **alpha;
  gt_assert(L);
  alpha = lua_newuserdata(L, sizeof *alpha);
  gt_assert(alpha);
  *alpha = gt_alphabet_new_protein();
  gt_assert(*alpha);
  luaL_getmetatable(L, ALPHABET_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int alphabet_lua_new_empty(lua_State *L)
{
  GtAlphabet **alpha;
  gt_assert(L);
  alpha = lua_newuserdata(L, sizeof *alpha);
  gt_assert(alpha);
  *alpha = gt_alphabet_new_empty();
  gt_assert(*alpha);
  luaL_getmetatable(L, ALPHABET_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int alphabet_lua_new_dna(lua_State *L)
{
  GtAlphabet **alpha;
  gt_assert(L);
  alpha = lua_newuserdata(L, sizeof *alpha);
  gt_assert(alpha);
  *alpha = gt_alphabet_new_dna();
  gt_assert(*alpha);
  luaL_getmetatable(L, ALPHABET_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int alphabet_lua_add_mapping(lua_State *L)
{
  GtAlphabet **alpha;
  const char *characters;
  alpha = check_alphabet(L, 1);
  characters = luaL_checkstring(L, 2);
  luaL_argcheck(L, strlen(characters), 2,
                "mapping must contain at least one character");
  gt_alphabet_add_mapping(*alpha, characters);
  return 0;
}

static int alphabet_lua_add_wildcard(lua_State *L)
{
  GtAlphabet **alpha;
  const char *wildcard;
  alpha = check_alphabet(L, 1);
  wildcard = luaL_checkstring(L, 2);
  luaL_argcheck(L, strlen(wildcard) == 1, 2,
                "wildcard string must have length 0");
  gt_alphabet_add_wildcard(*alpha, wildcard[0]);
  return 0;
}

static int alphabet_lua_decode(lua_State *L)
{
  GtAlphabet **alpha;
  unsigned int code;
  char character;
  alpha = check_alphabet(L, 1);
  code = luaL_checkinteger(L, 2);
  /* XXX: too restrictive, does not consider wildcards */
  luaL_argcheck(L, code < gt_alphabet_size(*alpha), 2, "invalid code");
  character = gt_alphabet_decode(*alpha, code);
  lua_pushlstring(L, &character, 1);
  return 1;
}

static int alphabet_lua_size(lua_State *L)
{
  GtAlphabet **alpha;
  unsigned int size;
  alpha = check_alphabet(L, 1);
  size = gt_alphabet_size(*alpha);
  lua_pushinteger(L, size);
  return 1;
}

static int alphabet_lua_delete(lua_State *L)
{
  GtAlphabet **alpha;
  alpha = check_alphabet(L, 1);
  gt_alphabet_delete(*alpha);
  return 0;
}

static const struct luaL_Reg alphabet_lib_f [] = {
  { "alphabet_new_dna", alphabet_lua_new_dna },
  { "alphabet_new_protein", alphabet_lua_new_protein },
  { "alphabet_new_empty", alphabet_lua_new_empty },
  { NULL, NULL }
};

static const struct luaL_Reg alphabet_lib_m [] = {
  { "add_mapping", alphabet_lua_add_mapping },
  { "add_wildcard", alphabet_lua_add_wildcard },
  { "decode", alphabet_lua_decode },
  { "size", alphabet_lua_size },
  { NULL, NULL }
};

int gt_lua_open_alphabet(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_newmetatable(L, ALPHABET_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, alphabet_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, alphabet_lib_m);
  lua_pop(L, 1);
  luaL_register(L, "gt", alphabet_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}

void gt_lua_alphabet_push(lua_State *L, GtAlphabet *alpha)
{
  GtAlphabet **alphaptr;
  gt_assert(L && alpha);
  alphaptr = lua_newuserdata(L, sizeof (*alphaptr));
  *alphaptr = alpha;
  luaL_getmetatable(L, ALPHABET_METATABLE);
  lua_setmetatable(L, -2);
}
