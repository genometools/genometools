/*
  Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2010 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "core/assert_api.h"
#include "core/ensure.h"
#include "core/unused_api.h"
#include "extended/luaserialize.h"

static int format_scalar(lua_State *L, GtStr *out, int index, bool table_key,
                         GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(!lua_istable(L, index));
  if (lua_isboolean(L, index))
  {
    int val;
    val = lua_toboolean(L, index);
    if (val)
      gt_str_append_cstr(out, "true");
    else
      gt_str_append_cstr(out, "false");
  }
  else if (lua_isnumber(L, index))
  {
    double val;
    val = lua_tonumber(L, index);
    gt_str_append_double(out, val, 10);
  }
  else if (lua_isstring(L, index))
  {
    const char *str;
    str = lua_tostring(L, index);
    /* use Lua's own "string.format() function to escape string */
    lua_getglobal(L, "string");
    lua_pushliteral(L, "format");
    lua_gettable(L, -2);
    lua_pushstring(L, "%q");
    lua_pushstring(L, str);
    lua_call(L, 2 ,1);
    str = lua_tostring(L, -1);
    lua_pop(L,2);
    /* table keys must be enclosed in square brackets */
    if (table_key)
      gt_str_append_cstr(out, "[");
    gt_str_append_cstr(out, str);
    if (table_key)
      gt_str_append_cstr(out, "]");
  }
  else {
    lua_pop(L, 2);
    gt_error_set(err, "expected boolean, number, or string");
    had_err = -1;
  }
  return had_err;
}

static int parse_table(lua_State *L, GtStr *out, int index, int level,
                       GtError *err)
{
  int GT_UNUSED rval, had_err = 0;
  gt_error_check(err);
  gt_assert(lua_istable(L, index));
  lua_pushnil(L);
  if (index < 0)
    index--;
  while (!had_err && (lua_next(L, index) != 0))
  {
    int i;
    for (i=0;i<level;i++)
      gt_str_append_cstr(out, "  ");
    rval = format_scalar(L, out, -2, true, NULL);
    gt_assert(!rval); /* cannot happen */
    gt_str_append_cstr(out, " = ");
    if (lua_istable(L, -1))
    {
      gt_str_append_cstr(out, "{\n");
      had_err = parse_table(L, out, -1, level+1, err);
      for (i=0;i<level;i++)
        gt_str_append_cstr(out, "  ");
      gt_str_append_cstr(out, "},\n");
    }
    else
    {
      had_err = format_scalar(L, out, -1, false, err);
      gt_str_append_cstr(out, ",\n");
    }
    lua_pop(L, 1);
  }
  return had_err;
}

int gt_lua_table_to_str(lua_State *L, GtStr *out, int index, GtError *err)
{
  int had_err;
#ifndef NDEBUG
  int stack_size = lua_gettop(L);
#endif
  gt_error_check(err);
  gt_assert(L && out && lua_istable(L, index));
  had_err = parse_table(L, out, index, 1, err);
  gt_assert(lua_gettop(L) == stack_size); /* make sure the stack doesn't grow */
  return had_err;
}

int gt_lua_serializer_unit_test(GtError *err)
{
  int had_err = 0;
  lua_State *L;
  GtStr *outstr  = gt_str_new();
  const char testtable[] = "config =\n"
  "{\n"
  "  gene = {\n"
  "    -- Color definitions\n"
  "    stroke             = {red=0.0, green=0.0, blue=0.0},\n"
  "    stroke_marked      = {red=1.0, green=0.0, blue=0.0},\n"
  "    fill               = {red=0.9, green=0.9, blue=1.0},\n"
  "    style              = \"box\",\n"
  "    -- Collapsing options\n"
  "    collapse_to_parent = false,\n"
  "    split_lines        = true,\n"
  "    -- Caption options\n"
  "    max_capt_show_width= nil,\n"
  "    -- Display this track only if the viewport is not wider than this\n"
  "    -- number of nucleotides. Set to 0 to disable type track.\n"
  "    max_show_width     = nil,\n"
  "    -- Limit the number of tracks\n"
  "    --    max_num_lines      = 10,\n"
  "  },\n"
  "}\n"
  "";
  gt_error_check(err);

  L = luaL_newstate();
  gt_ensure(had_err, L);
  luaL_openlibs(L);
  had_err = luaL_loadbuffer(L, testtable, sizeof (testtable)-1, "t") ||
              lua_pcall(L, 0, 0, 0);
  if (!had_err)
  {
    lua_getglobal(L, "config");
    gt_str_append_cstr(outstr, "config = {\n");
    gt_lua_table_to_str(L, outstr, -1, err);
    gt_str_append_cstr(outstr, "}");
    had_err = luaL_loadbuffer(L, gt_str_get(outstr),
                              gt_str_length(outstr), "t2") ||
              lua_pcall(L, 0, 0, 0);
    if (!had_err)
    {
      gt_ensure(had_err, gt_str_length(outstr) > 0);
      lua_getglobal(L, "config");
      gt_ensure(had_err, lua_istable(L, -1));
      lua_getfield(L, -1, "gene");
      gt_ensure(had_err, lua_istable(L, -1));
      lua_getfield(L, -1, "stroke_marked");
      gt_ensure(had_err, lua_istable(L, -1));
      lua_getfield(L, -1, "red");
      gt_ensure(had_err, lua_isnumber(L, -1));
      gt_ensure(had_err, 1.0 == lua_tonumber(L, -1));
      lua_pop(L, 1);
      lua_getfield(L, -1, "blue");
      gt_ensure(had_err, lua_isnumber(L, -1));
      gt_ensure(had_err, 0.0 == lua_tonumber(L, -1));
      lua_pop(L, 2);
      lua_getfield(L, -1, "style");
      gt_ensure(had_err, lua_isstring(L, -1));
      gt_ensure(had_err, strcmp("box",lua_tostring(L, -1)) == 0);
      lua_pop(L, 1);
      lua_getfield(L, -1, "collapse_to_parent");
      gt_ensure(had_err, lua_isboolean(L, -1));
      gt_ensure(had_err, lua_toboolean(L, -1) == false);
      lua_pop(L, 1);
      lua_getfield(L, -1, "max_show_width");
      gt_ensure(had_err, lua_isnil(L, -1));
      lua_pop(L, 3);
    }
  }
  gt_str_delete(outstr);
  lua_close(L);
  return had_err;
}
