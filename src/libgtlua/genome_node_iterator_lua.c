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

#include "lauxlib.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtlua/genome_node_lua.h"
#include "libgtlua/genome_node_iterator_lua.h"

#define GENOME_NODE_ITERATOR_METATABLE  "GenomeTools_genome.node_iterator"
#define check_genome_node_iterator(L, POS) \
        (GenomeNodeIterator**) \
        luaL_checkudata(L, POS, GENOME_NODE_ITERATOR_METATABLE)

static int genome_node_iterator_lua_new(lua_State *L)
{
  GenomeNodeIterator **gni;
  GenomeNode **gn;
  assert(L);
  gn = check_genome_node(L, 1);
  gni = lua_newuserdata(L, sizeof (GenomeNodeIterator*));
  assert(gni);
  *gni = genome_node_iterator_new(*gn);
  luaL_getmetatable(L, GENOME_NODE_ITERATOR_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int genome_node_iterator_lua_new_direct(lua_State *L)
{
  GenomeNodeIterator **gni;
  GenomeNode **gn;
  assert(L);
  gn = check_genome_node(L, 1);
  gni = lua_newuserdata(L, sizeof (GenomeNodeIterator*));
  assert(gni);
  *gni = genome_node_iterator_new_direct(*gn);
  luaL_getmetatable(L, GENOME_NODE_ITERATOR_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int genome_node_iterator_lua_next(lua_State *L)
{
  GenomeNodeIterator **gni;
  GenomeNode *gn;
  gni = check_genome_node_iterator(L, 1);
  gn = genome_node_iterator_next(*gni);
  if (gn)
    genome_node_lua_push(L, genome_node_rec_ref(gn));
  else
    lua_pushnil(L);
  return 1;
}

static int genome_node_iterator_lua_delete(lua_State *L)
{
  GenomeNodeIterator **genome_node_iterator;
  genome_node_iterator = check_genome_node_iterator(L, 1);
  genome_node_iterator_delete(*genome_node_iterator);
  return 0;
}

static const struct luaL_Reg genome_node_iterator_lib_f [] = {
  { "genome_node_iterator_new", genome_node_iterator_lua_new },
  { "genome_node_iterator_new_direct", genome_node_iterator_lua_new_direct },
  { NULL, NULL }
};

static const struct luaL_Reg genome_node_iterator_lib_m [] = {
  { "next", genome_node_iterator_lua_next },
  { NULL, NULL }
};

int luaopen_genome_node_iterator(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, GENOME_NODE_ITERATOR_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, genome_node_iterator_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, genome_node_iterator_lib_m);
  luaL_register(L, "gt", genome_node_iterator_lib_f);
  return 1;
}
