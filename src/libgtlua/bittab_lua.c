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
#include "libgtcore/bittab.h"
#include "libgtcore/error.h"
#include "libgtext/luahelper.h"
#include "libgtlua/bittab_lua.h"
#include "libgtlua/gtcore_lua.h"

#define BITTAB_METATABLE  "GenomeTools.bittab"
#define check_bittab(L, POS) \
        (Bittab**) luaL_checkudata(L, POS, BITTAB_METATABLE)

static int bittab_lua_new(lua_State *L)
{
  long num_of_bits;
  Bittab **bittab;
  assert(L);
  num_of_bits = luaL_checklong(L, 1);
  luaL_argcheck(L, num_of_bits > 0, 1, "must be > 0");
  bittab = lua_newuserdata(L, sizeof (Bittab*));
  assert(bittab);
  *bittab = bittab_new(num_of_bits);
  luaL_getmetatable(L, BITTAB_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static void get_bittab_and_bit(lua_State *L, Bittab ***bittab, long *bit)
{
  *bittab = check_bittab(L, 1);
  *bit = luaL_checklong(L, 2);
  luaL_argcheck(L, *bit >= 0, 2, "bit number too small");
  luaL_argcheck(L, *bit < bittab_size(**bittab), 2, "bit number too large");
}

static void get_two_bittabs(lua_State *L, Bittab ***bt1, Bittab ***bt2)
{
 *bt1 = check_bittab(L, 1);
 *bt2 = check_bittab(L, 2);
 luaL_argcheck(L, bittab_size(**bt1) == bittab_size(**bt2), 1, "bittabs have "
               "different sizes");
}

static int bittab_lua_set_bit(lua_State *L)
{
  Bittab **bittab;
  long bit;
  get_bittab_and_bit(L, &bittab, &bit);
  bittab_set_bit(*bittab, bit);
  return 0;
}

static int bittab_lua_unset_bit(lua_State *L)
{
  Bittab **bittab;
  long bit;
  get_bittab_and_bit(L, &bittab, &bit);
  bittab_unset_bit(*bittab, bit);
  return 0;
}

static int bittab_lua_complement(lua_State *L)
{
  Bittab **dest, **src;
  get_two_bittabs(L, &dest, &src);
  bittab_complement(*dest, *src);
  return 0;
}

static int bittab_lua_equal(lua_State *L)
{
  Bittab **dest, **src;
  get_two_bittabs(L, &dest, &src);
  bittab_equal(*dest, *src);
  return 0;
}

static int bittab_lua_and_equal(lua_State *L)
{
  Bittab **dest, **src;
  get_two_bittabs(L, &dest, &src);
  bittab_and_equal(*dest, *src);
  return 0;
}

static int bittab_lua_bit_is_set(lua_State *L)
{
  Bittab **bittab;
  long bit;
  get_bittab_and_bit(L, &bittab, &bit);
  lua_pushboolean(L, bittab_bit_is_set(*bittab, bit));
  return 1;
}

static int bittab_lua_delete(lua_State *L)
{
  Bittab **bittab;
  bittab = check_bittab(L, 1);
  bittab_delete(*bittab);
  return 0;
}

static const struct luaL_Reg bittab_lib_f [] = {
  { "bittab_new", bittab_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg bittab_lib_m [] = {
  { "set_bit", bittab_lua_set_bit },
  { "unset_bit", bittab_lua_unset_bit },
  { "complement", bittab_lua_complement },
  { "equal", bittab_lua_equal},
  { "and_equal", bittab_lua_and_equal},
  { "bit_is_set", bittab_lua_bit_is_set},
  { NULL, NULL }
};

int luaopen_bittab(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, BITTAB_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, bittab_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, bittab_lib_m);
  luaL_register(L, "gt", bittab_lib_f);
  return 1;
}
