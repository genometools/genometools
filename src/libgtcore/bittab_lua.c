/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "lauxlib.h"
#include "gtlua.h"
#include "libgtcore/bittab.h"
#include "libgtcore/bittab_lua.h"
#include "libgtcore/env.h"
#include "libgtcore/gtcore_lua.h"

#define BITTAB_METATABLE  "GenomeTools.bittab"
#define checkbittab(L, POS) \
        (Bittab**) luaL_checkudata(L, POS, BITTAB_METATABLE);

static int bittab_lua_new(lua_State *L)
{
  long num_of_bits;
  Bittab **bittab;
  Env *env;
  assert(L);
  num_of_bits = luaL_checklong(L, 1);
  luaL_argcheck(L, num_of_bits > 0, 1, "must be > 0");
  env = get_env_from_registry(L);
  bittab = lua_newuserdata(L, sizeof (Bittab**));
  assert(bittab);
  *bittab = bittab_new(num_of_bits, env);
  luaL_getmetatable(L, BITTAB_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static void get_bittab_and_bit(lua_State *L, Bittab ***bittab, long *bit)
{
  *bittab = checkbittab(L, 1);
  *bit = luaL_checklong(L, 2);
  luaL_argcheck(L, *bit < bittab_size(**bittab), 2, "bit number too large");
}

static void get_two_bittabs(lua_State *L, Bittab ***bt1, Bittab ***bt2)
{
 *bt1 = checkbittab(L, 1);
 *bt2 = checkbittab(L, 2);
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

static int bittab_lua_delete(lua_State *L)
{
  Bittab **bittab;
  Env *env = get_env_from_registry(L);
  bittab = checkbittab(L, 1);
  bittab_delete(*bittab, env);
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
