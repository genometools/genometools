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
#include "libgtext/gff3_visitor.h"
#include "libgtext/luahelper.h"
#include "libgtlua/genome_visitor_lua.h"

static int gff3_visitor_lua_new(lua_State *L)
{
  GenomeVisitor **gv;
  assert(L);
  /* construct object */
  gv = lua_newuserdata(L, sizeof (GenomeVisitor*));
  *gv = gff3_visitor_new(NULL);
  assert(*gv);
  luaL_getmetatable(L, GENOME_VISITOR_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int genome_visitor_lua_delete(lua_State *L)
{
  GenomeVisitor **gv = check_genome_visitor(L, 1);
  genome_visitor_delete(*gv);
  return 0;
}

static const struct luaL_Reg genome_visitor_lib_f [] = {
  { "gff3_visitor_new", gff3_visitor_lua_new },
  { NULL, NULL }
};

int luaopen_genome_visitor(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, GENOME_VISITOR_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, genome_visitor_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, "gt", genome_visitor_lib_f);
  return 1;
}
