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

#include <string.h>
#include "lauxlib.h"
#include "core/assert_api.h"
#include "core/cstr_api.h"
#include "core/gtdatapath.h"
#include "core/ma.h"
#include "extended/luahelper.h"

/* key used to store the Env object in the Lua registry */
#define ENV_KEY env_new

int gt_lua_set_modules_path(lua_State *L, GtError *err)
{
  GtStr *modules_path = NULL, *external_modules_path = NULL,
         *package_path = NULL;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(L);
  if (!(modules_path = gt_get_gtdata_path(gt_error_get_progname(err), err)))
    had_err = -1;
  if (!had_err) {
    external_modules_path = gt_str_clone(modules_path);
    gt_str_append_cstr(modules_path, "/modules/?.lua");
    gt_str_append_cstr(external_modules_path, "/modules/external/?.lua");
    lua_getglobal(L, "package");
    gt_assert(lua_istable(L, -1));
    lua_getfield(L, -1, "path");
    gt_assert(lua_isstring(L, -1));
    package_path = gt_str_new_cstr(lua_tostring(L, -1));
    lua_pop(L, 1);
    gt_str_append_char(package_path, ';');
    gt_str_append_str(package_path, modules_path);
    gt_str_append_char(package_path, ';');
    gt_str_append_str(package_path, external_modules_path);
    lua_pushstring(L, gt_str_get(package_path));
    lua_setfield(L, -2, "path");
    lua_pop(L, 1);
  }
  gt_str_delete(package_path);
  gt_str_delete(modules_path);
  gt_str_delete(external_modules_path);
  return had_err;
}

void gt_lua_set_arg(lua_State *L, const char *argv_0, const char **argv)
{
  lua_Integer n = 0;
  gt_assert(L && argv_0);
  /* create table */
  lua_newtable(L);
  /* set arg[0] */
  lua_pushinteger(L, 0);
  lua_pushstring(L, argv_0);
  lua_rawset(L, -3);
  /* set other arguments */
  while (argv[n]) {
    lua_pushinteger(L, n+1);
    lua_pushstring(L, argv[n]);
    lua_rawset(L, -3);
    n++;
  }
  /* register table globally */
  lua_setglobal(L, "arg");
}

void gt_lua_export_metatable(lua_State *L, const char *metatable_desc)
{
  char *dot, *mt;
  gt_assert(L && metatable_desc);
  mt = gt_cstr_dup(metatable_desc);
  dot = strchr(mt, '.');
  gt_assert(dot);
  *dot = '_';
  lua_setglobal(L, mt);
  gt_free(mt);
}

void gt_lua_push_strarray_as_table(lua_State *L, GtStrArray *sa)
{
  GtUword i;
  gt_assert(L && sa);
  lua_newtable(L);
  for (i = 0; i < gt_str_array_size(sa); i++) {
    lua_pushinteger(L, i+1); /* in Lua we index from 1 on */
    lua_pushstring(L, gt_str_array_get(sa, i));
    lua_rawset(L, -3);
  }
}

int gt_lua_get_table_as_strarray(lua_State *L, int index, GtStrArray *outarray,
                                 GtError *err)
{
  int had_err = 0;
  gt_assert(lua_istable(L, index));
  lua_pushnil(L);
  while (!had_err && (lua_next(L, index) != 0))
  {
    if (!lua_isstring(L, -1)) {
      had_err = -1;
      gt_error_set(err, "table contains non-string value!");
      break;
    }
    gt_str_array_add_cstr(outarray, lua_tostring(L, -1));
    lua_pop(L, 1);
  }
  return had_err;
}

int gt_lua_error(lua_State *L, GtError *err)
{
  gt_assert(L && err);
  gt_assert(gt_error_is_set(err));
  luaL_where(L, 1);
  lua_pushstring(L, gt_error_get(err));
  gt_error_delete(err);
  lua_concat(L, 2);
  return lua_error(L);
}

void* gt_lua_try_checkudata(lua_State *L, int ud, const char *tname)
{
  void *p;
  gt_assert(L && tname);
  p = lua_touserdata(L, ud);
  if (p != NULL) {  /* value is a userdata? */
    if (lua_getmetatable(L, ud)) {  /* does it have a metatable? */
      lua_getfield(L, LUA_REGISTRYINDEX, tname);  /* get correct metatable */
      if (lua_rawequal(L, -1, -2)) {  /* does it have the correct mt? */
        lua_pop(L, 2);  /* remove both metatables */
        return p;
      }
    }
  }
  return NULL;  /* to avoid warnings */
}
