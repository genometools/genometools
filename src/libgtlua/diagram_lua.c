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

#ifdef LIBGTVIEW

#include "lauxlib.h"
#include "libgtext/luahelper.h"
#include "libgtlua/diagram_lua.h"
#include "libgtlua/feature_index_lua.h"
#include "libgtlua/genome_node_lua.h"
#include "libgtlua/range_lua.h"
#include "libgtview/feature_index.h"
#include "libgtview/diagram.h"
#include "libgtview/luaconfig.h"

static int diagram_lua_new(lua_State *L)
{
  Diagram **diagram;
  FeatureIndex **feature_index;
  Range *range;
  const char *seqid;
  Config *config;
  /* get feature index */
  feature_index = check_feature_index(L, 1);
  /* get seqid */
  seqid       = luaL_checkstring(L, 2);
  luaL_argcheck(L, feature_index_has_seqid(*feature_index, seqid),
                2, "feature index does not contain the given sequence id");
  /* get range */
  range = check_range(L, 3);
  /* create diagram */
  config = lua_get_config_from_registry(L);
  diagram = lua_newuserdata(L, sizeof (Diagram*));
  assert(diagram);
  *diagram = diagram_new(*feature_index, seqid, range, config);
  luaL_getmetatable(L, DIAGRAM_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int diagram_lua_delete(lua_State *L)
{
  Diagram **diagram;
  diagram = check_diagram(L, 1);
  diagram_delete(*diagram);
  return 0;
}

static const struct luaL_Reg diagram_lib_f [] = {
  { "diagram_new", diagram_lua_new },
  { NULL, NULL }
};

int luaopen_diagram(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, DIAGRAM_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, diagram_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, "gt", diagram_lib_f);
  return 1;
}

#endif
