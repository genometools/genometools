/*
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
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

#include "lauxlib.h"
#include "extended/feature_node_iterator_api.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/feature_node_iterator_lua.h"

#define GENOME_NODE_ITERATOR_METATABLE  "GenomeTools.feature_node_iterator"
#define check_gt_feature_node_iterator(L, POS) \
        (GtFeatureNodeIterator**) \
        luaL_checkudata(L, POS, GENOME_NODE_ITERATOR_METATABLE)

static int feature_node_iterator_lua_new(lua_State *L)
{
  GtFeatureNodeIterator **fni;
  GtFeatureNode **fn;
  gt_assert(L);
  fn = (GtFeatureNode**) check_genome_node(L, 1);
  fni = lua_newuserdata(L, sizeof (GtFeatureNodeIterator*));
  gt_assert(fni);
  *fni = gt_feature_node_iterator_new(*fn);
  luaL_getmetatable(L, GENOME_NODE_ITERATOR_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int feature_node_iterator_lua_new_direct(lua_State *L)
{
  GtFeatureNodeIterator **fni;
  GtFeatureNode **fn;
  gt_assert(L);
  fn = (GtFeatureNode**) check_genome_node(L, 1);
  fni = lua_newuserdata(L, sizeof (GtFeatureNodeIterator*));
  gt_assert(fni);
  *fni = gt_feature_node_iterator_new_direct(*fn);
  luaL_getmetatable(L, GENOME_NODE_ITERATOR_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int feature_node_iterator_lua_next(lua_State *L)
{
  GtFeatureNodeIterator **fni;
  GtFeatureNode *fn;
  fni = check_gt_feature_node_iterator(L, 1);
  fn = gt_feature_node_iterator_next(*fni);
  if (fn)
    gt_lua_genome_node_push(L, gt_genome_node_ref((GtGenomeNode*) fn));
  else
    lua_pushnil(L);
  return 1;
}

static int feature_node_iterator_lua_delete(lua_State *L)
{
  GtFeatureNodeIterator **gt_feature_node_iterator;
  gt_feature_node_iterator = check_gt_feature_node_iterator(L, 1);
  gt_feature_node_iterator_delete(*gt_feature_node_iterator);
  return 0;
}

static const struct luaL_Reg feature_node_iterator_lib_f [] = {
  { "feature_node_iterator_new", feature_node_iterator_lua_new },
  { "feature_node_iterator_new_direct", feature_node_iterator_lua_new_direct },
  { NULL, NULL }
};

static const struct luaL_Reg feature_node_iterator_lib_m [] = {
  { "next", feature_node_iterator_lua_next },
  { NULL, NULL }
};

int gt_lua_open_feature_node_iterator(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_newmetatable(L, GENOME_NODE_ITERATOR_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, feature_node_iterator_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, feature_node_iterator_lib_m);
  lua_pop(L, 1);
  luaL_register(L, "gt", feature_node_iterator_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}
