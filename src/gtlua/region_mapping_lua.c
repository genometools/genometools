/*
  Copyright (c) 2008 Gordon Gremme <gordon@gremme.org>
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
#include "extended/luahelper.h"
#include "extended/region_mapping.h"
#include "gtlua/region_mapping_lua.h"
#include "gtlua/gtcore_lua.h"

static int region_mapping_lua_new_seqfile_matchdesc(lua_State *L)
{
  const char *seqfilename;
  GtStrArray *seqfile;
  GtRegionMapping **region_mapping;
  gt_assert(L);
  seqfilename = luaL_checkstring(L, 1);
  region_mapping = lua_newuserdata(L, sizeof (GtRegionMapping*));
  gt_assert(region_mapping);
  seqfile = gt_str_array_new();
  gt_str_array_add_cstr(seqfile, seqfilename);
  /* XXX: make second and third parameter available */
  *region_mapping = gt_region_mapping_new_seqfiles(seqfile, true, false);
  gt_str_array_delete(seqfile);
  luaL_getmetatable(L, REGION_MAPPING_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int region_mapping_lua_new_seqfile_matchdescstart(lua_State *L)
{
  const char *seqfilename;
  GtStrArray *seqfile;
  GtRegionMapping **region_mapping;
  gt_assert(L);
  seqfilename = luaL_checkstring(L, 1);
  region_mapping = lua_newuserdata(L, sizeof (GtRegionMapping*));
  gt_assert(region_mapping);
  seqfile = gt_str_array_new();
  gt_str_array_add_cstr(seqfile, seqfilename);
  /* XXX: make second and third parameter available */
  *region_mapping = gt_region_mapping_new_seqfiles(seqfile, true, false);
  gt_region_mapping_enable_match_desc_start(*region_mapping);
  gt_str_array_delete(seqfile);
  luaL_getmetatable(L, REGION_MAPPING_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int region_mapping_lua_new_seqfile_usedesc(lua_State *L)
{
  const char *seqfilename;
  GtStrArray *seqfile;
  GtRegionMapping **region_mapping;
  gt_assert(L);
  seqfilename = luaL_checkstring(L, 1);
  region_mapping = lua_newuserdata(L, sizeof (GtRegionMapping*));
  gt_assert(region_mapping);
  seqfile = gt_str_array_new();
  gt_str_array_add_cstr(seqfile, seqfilename);
  /* XXX: make second and third parameter available */
  *region_mapping = gt_region_mapping_new_seqfiles(seqfile, false, true);
  gt_str_array_delete(seqfile);
  luaL_getmetatable(L, REGION_MAPPING_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int region_mapping_lua_new_seqfile(lua_State *L)
{
  return region_mapping_lua_new_seqfile_matchdesc(L);
}

static int region_mapping_lua_delete(lua_State *L)
{
  GtRegionMapping **region_mapping;
  region_mapping = check_region_mapping(L, 1);
  gt_region_mapping_delete(*region_mapping);
  return 0;
}

static const struct luaL_Reg region_mapping_lib_f [] = {
  { "region_mapping_new_seqfile", region_mapping_lua_new_seqfile },
  { "region_mapping_new_seqfile_matchdesc",
                                     region_mapping_lua_new_seqfile_matchdesc },
  { "region_mapping_new_seqfile_matchdescstart",
                                region_mapping_lua_new_seqfile_matchdescstart },
  { "region_mapping_new_seqfile_usedesc",
                                       region_mapping_lua_new_seqfile_usedesc },
  { NULL, NULL }
};

int gt_lua_open_region_mapping(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_newmetatable(L, REGION_MAPPING_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, region_mapping_lua_delete);
  lua_settable(L, -3);
  lua_pop(L, 1);
  /* register functions */
  luaL_register(L, "gt", region_mapping_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}
