/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
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

#include <assert.h>
#include <string.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "libgtcore/ensure.h"
#include "libgtcore/unused.h"
#include "libgtview/luaserialize.h"

static int format_scalar(lua_State *L, Str *out, int index, bool table_key)
{
  int had_err = 0;
  assert(!lua_istable(L ,index));
  if (lua_isboolean(L, index))
  {
    int val;
    val = lua_toboolean(L, index);
    if (val)
      str_append_cstr(out, "true");
    else
      str_append_cstr(out, "false");
  }
  else if (lua_isnumber(L, index))
  {
    double val;
    val = lua_tonumber(L, index);
    str_append_double(out, val, 10);
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
      str_append_cstr(out, "[");
    str_append_cstr(out, str);
    if (table_key)
      str_append_cstr(out, "]");
  } else had_err = -1;
  return had_err;
}

static int parse_table(lua_State *L, Str *out, int index, int level)
{
  int had_err = 0;
  if (!lua_istable(L, index))
    return -1;
  lua_pushnil(L);
  if (index < 0)
    index--;
  while (!had_err && (lua_next(L, index) != 0))
  {
    int i;
    for (i=0;i<level;i++)
      str_append_cstr(out, "  ");
    format_scalar(L, out, -2, true);
    str_append_cstr(out, " = ");
    if (lua_istable(L, -1))
    {
      str_append_cstr(out, "{\n");
      had_err = parse_table(L, out, -1, level+1);
      for (i=0;i<level;i++)
        str_append_cstr(out, "  ");
      str_append_cstr(out, "},\n");
    }
    else
    {
      had_err = format_scalar(L, out, -1, false);
      str_append_cstr(out, ",\n");
    }
    lua_pop(L, 1);
  }
  return had_err;
}

int lua_table_to_str(lua_State *L, Str *out, int index)
{
  int had_err = 0;
  assert(L && out);
  had_err = parse_table(L, out, index, 1);
  return had_err;
}

int lua_serializer_unit_test(Error *err)
{
  int had_err = 0;
  lua_State *L;
  Str *outstr  = str_new();
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
  error_check(err);

  L = luaL_newstate();
  ensure(had_err, L);
  luaL_openlibs(L);
  had_err = luaL_loadbuffer(L, testtable, sizeof (testtable)-1, "t") ||
              lua_pcall(L, 0, 0, 0);
  if (!had_err)
  {
    lua_getglobal(L, "config");
    str_append_cstr(outstr, "config = {\n");
    lua_table_to_str(L, outstr, -1);
    str_append_cstr(outstr, "}");
    had_err = luaL_loadbuffer(L, str_get(outstr), str_length(outstr), "t2") ||
                lua_pcall(L, 0, 0, 0);
    if (!had_err)
    {
      ensure(had_err, str_length(outstr) > 0);
      lua_getglobal(L, "config");
      ensure(had_err, lua_istable(L, -1));
      lua_getfield(L, -1, "gene");
      ensure(had_err, lua_istable(L, -1));
      lua_getfield(L, -1, "stroke_marked");
      ensure(had_err, lua_istable(L, -1));
      lua_getfield(L, -1, "red");
      ensure(had_err, lua_isnumber(L, -1));
      ensure(had_err, 1.0 == lua_tonumber(L, -1));
      lua_pop(L, 1);
      lua_getfield(L, -1, "blue");
      ensure(had_err, lua_isnumber(L, -1));
      ensure(had_err, 0.0 == lua_tonumber(L, -1));
      lua_pop(L, 2);
      lua_getfield(L, -1, "style");
      ensure(had_err, lua_isstring(L, -1));
      ensure(had_err, strcmp("box",lua_tostring(L, -1)) == 0);
      lua_pop(L, 1);
      lua_getfield(L, -1, "collapse_to_parent");
      ensure(had_err, lua_isboolean(L, -1));
      ensure(had_err, lua_toboolean(L, -1) == false);
      lua_pop(L, 1);
      lua_getfield(L, -1, "max_show_width");
      ensure(had_err, lua_isnil(L, -1));
      lua_pop(L, 3);
    }
  }
  str_delete(outstr);
  lua_close(L);
  return had_err;
}
