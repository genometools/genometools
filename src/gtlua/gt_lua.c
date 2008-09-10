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

#include <assert.h>
#include "gtlua/gtcore_lua.h"
#include "gtlua/gtext_lua.h"
#include "gtlua/gt_lua.h"

#ifndef WITHOUT_CAIRO
#include "gtlua/annotationsketch_lua.h"
#endif

/* key used to store the GT_TypeFactory object in the Lua registry */
#define FEATURE_TYPE_FACTORY_KEY gt_type_factory_create_gft

void gt_lua_put_type_factory_in_registry(lua_State *L,
                                         GT_TypeFactory *feature_type_factory)
{
  assert(L && feature_type_factory);
  lua_pushlightuserdata(L, FEATURE_TYPE_FACTORY_KEY); /* push the key */
  lua_pushlightuserdata(L, feature_type_factory); /* push the value */
  lua_rawset(L, LUA_REGISTRYINDEX); /* store feature type factory in registry */
}

GT_TypeFactory* gt_lua_get_type_factory_from_registry(lua_State *L)
{
  GT_TypeFactory *feature_type_factory;
  assert(L);
  lua_pushlightuserdata(L, FEATURE_TYPE_FACTORY_KEY);
  lua_rawget(L, LUA_REGISTRYINDEX);
  assert(lua_islightuserdata(L, -1));
  feature_type_factory = lua_touserdata(L, -1);
  lua_pop(L, 1);
  return feature_type_factory;
}

int gt_lua_openlib(lua_State *L)
{
  assert(L);
  gt_lua_open_core(L); /* open core */
  gt_lua_open_extended(L); /* open extended */
#ifndef WITHOUT_CAIRO
  gt_lua_open_annotationsketch(L); /* open annotationsketch */
#endif
  return 1;
}
