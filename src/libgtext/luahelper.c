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

#include <assert.h>
#include <string.h>
#include "lauxlib.h"
#include "libgtcore/cstr.h"
#include "libgtcore/gtdatapath.h"
#include "libgtcore/ma.h"
#include "libgtext/luahelper.h"

/* key used to store the Env object in the Lua registry */
#define ENV_KEY env_new

int lua_set_modules_path(lua_State *L, Error *err)
{
  Str *modules_path = NULL, *external_modules_path = NULL, *package_path = NULL;
  int had_err = 0;
  error_check(err);
  assert(L);
  if (!(modules_path = gtdata_get_path(error_get_progname(err), err)))
    had_err = -1;
  if (!had_err) {
    external_modules_path = str_clone(modules_path);
    str_append_cstr(modules_path, "/modules/?.lua");
    str_append_cstr(external_modules_path, "/modules/external/?.lua");
    lua_getglobal(L, "package");
    assert(lua_istable(L, -1));
    lua_getfield(L, -1, "path");
    assert(lua_isstring(L, -1));
    package_path = str_new_cstr(lua_tostring(L, -1));
    lua_pop(L, 1);
    str_append_char(package_path, ';');
    str_append_str(package_path, modules_path);
    str_append_char(package_path, ';');
    str_append_str(package_path, external_modules_path);
    lua_pushstring(L, str_get(package_path));
    lua_setfield(L, -2, "path");
    lua_pop(L, 1);
  }
  str_delete(package_path);
  str_delete(modules_path);
  str_delete(external_modules_path);
  return had_err;
}

void lua_set_arg(lua_State *L, const char *argv_0, const char **argv)
{
  lua_Integer n = 0;
  assert(L && argv_0);
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

void lua_export_metatable(lua_State *L, const char *metatable_desc)
{
  char *dot, *mt;
  assert(L && metatable_desc);
  mt = cstr_dup(metatable_desc);
  dot = strchr(mt, '.');
  assert(dot);
  *dot = '_';
  lua_setglobal(L, mt);
  ma_free(mt);
}

void lua_push_strarray_as_table(lua_State *L, StrArray *sa)
{
  unsigned long i;
  assert(L && sa);
  lua_newtable(L);
  for (i = 0; i < strarray_size(sa); i++) {
    lua_pushinteger(L, i+1); /* in Lua we index from 1 on */
    lua_pushstring(L, strarray_get(sa, i));
    lua_rawset(L, -3);
  }
}

int lua_gt_error(lua_State *L, Error *err)
{
  assert(L && err);
  assert(error_is_set(err));
  luaL_where(L, 1);
  lua_pushstring(L, error_get(err));
  error_delete(err);
  lua_concat(L, 2);
  return lua_error(L);
}
