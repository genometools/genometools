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
#include "libgtext/luahelper.h"
#include "libgtlua/genome_visitor_lua.h"
#include "libgtlua/feature_index_lua.h"
#include "libgtlua/feature_visitor_lua.h"
#include "libgtview/feature_visitor.h"

static int feature_visitor_lua_new(lua_State *L)
{
  GenomeVisitor **feature_visitor;
  FeatureIndex **feature_index;
  feature_visitor = lua_newuserdata(L, sizeof (GenomeVisitor*));
  assert(feature_visitor);
  feature_index = check_feature_index(L, 1);
  *feature_visitor = feature_visitor_new(*feature_index);
  luaL_getmetatable(L, GENOME_VISITOR_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static const struct luaL_Reg feature_visitor_lib_f [] = {
  { "feature_visitor_new", feature_visitor_lua_new },
  { NULL, NULL }
};

int luaopen_feature_visitor(lua_State *L)
{
  assert(L);
  luaL_register(L, "gt", feature_visitor_lib_f);
  return 1;
}
