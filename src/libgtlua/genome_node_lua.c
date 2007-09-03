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
#include "libgtext/genome_node.h"
#include "libgtlua/genome_node_lua.h"
#include "libgtlua/genome_visitor_lua.h"
#include "libgtlua/gt_lua.h"

#define check_genome_node(L) \
              (GenomeNode**) luaL_checkudata(L, 1, GENOME_NODE_METATABLE);

static int genome_feature_lua_new(lua_State *L)
{
  GenomeNode **gf;
  GenomeFeatureType type;
  Range range;
  Strand strand;
  const char *type_str, *strand_str;
  size_t length;
  Str *filename;
  Env *env;
  assert(L);
  /* get/check parameters */
  type_str = luaL_checkstring(L, 1);
  luaL_argcheck(L, !genome_feature_type_get(&type, type_str), 1,
                "invalid feature type");

  range.start = luaL_checklong(L, 2);
  range.end   = luaL_checklong(L, 3);
  luaL_argcheck(L, range.start <= range.end, 2, "must be <= endpos");

  strand_str = luaL_checklstring(L, 4, &length);
  luaL_argcheck(L, length == 1, 4, "strand string must have length 1");
  luaL_argcheck(L, (strand = strand_get(strand_str[0])) != NUM_OF_STRAND_TYPES,
                4, "invalid strand");
  /* construct object */
  env = get_env_from_registry(L);
  gf = lua_newuserdata(L, sizeof (GenomeNode**));
  filename = str_new_cstr("Lua", env);
  *gf = genome_feature_new(type, range, strand, filename, 0, env);
  str_delete(filename, env);
  assert(*gf);
  luaL_getmetatable(L, GENOME_NODE_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int genome_node_lua_get_filename(lua_State *L)
{
  GenomeNode **gn = check_genome_node(L);
  lua_pushstring(L, genome_node_get_filename(*gn));
  return 1;
}

static int genome_node_lua_accept(lua_State *L)
{
  GenomeNode **gn;
  GenomeVisitor **gv;
  Env *env;
  gn = check_genome_node(L);
  gv = check_genome_visitor(L, 2);
  env = get_env_from_registry(L);
  env_error_check(env);
  if (genome_node_accept(*gn, *gv, env))
    return luagt_error(L, env);
  return 0;
}

static int genome_node_lua_delete(lua_State *L)
{
  GenomeNode **gn = check_genome_node(L);
  Env *env;
  env = get_env_from_registry(L);
  genome_node_rec_delete(*gn, env);
  return 0;
}

static const struct luaL_Reg genome_node_lib_f [] = {
  { "genome_feature_new", genome_feature_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg genome_node_lib_m [] = {
  { "get_filename", genome_node_lua_get_filename },
  { "accept", genome_node_lua_accept },
  { NULL, NULL }
};

int luaopen_genome_node(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, GENOME_NODE_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, genome_node_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, genome_node_lib_m);
  luaL_register(L, "gt", genome_node_lib_f);
  return 1;
}

void genome_node_lua_push(lua_State *L, GenomeNode *gn)
{
  GenomeNode **gn_lua;
  assert(L && gn);
  gn_lua = lua_newuserdata(L, sizeof (GenomeNode**));
  *gn_lua = gn;
  luaL_getmetatable(L, GENOME_NODE_METATABLE);
  lua_setmetatable(L, -2);
}
