/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/region_mapping.h"
#include "libgtlua/region_mapping_lua.h"
#include "libgtlua/gtcore_lua.h"

static int region_mapping_lua_new_seqfile(lua_State *L)
{
  const char *seqfilename;
  Str *seqfile;
  RegionMapping **region_mapping;
  assert(L);
  seqfilename = luaL_checkstring(L, 1);
  region_mapping = lua_newuserdata(L, sizeof (RegionMapping*));
  assert(region_mapping);
  seqfile = str_new_cstr(seqfilename);
  *region_mapping = region_mapping_new_seqfile(seqfile);
  str_delete(seqfile);
  luaL_getmetatable(L, REGION_MAPPING_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int region_mapping_lua_delete(lua_State *L)
{
  RegionMapping **region_mapping;
  region_mapping = check_region_mapping(L, 1);
  region_mapping_delete(*region_mapping);
  return 0;
}

static const struct luaL_Reg region_mapping_lib_f [] = {
  { "region_mapping_new_seqfile", region_mapping_lua_new_seqfile },
  { NULL, NULL }
};

int luaopen_region_mapping(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, REGION_MAPPING_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, region_mapping_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, "gt", region_mapping_lib_f);
  return 1;
}
