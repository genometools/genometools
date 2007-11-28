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

#include <assert.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "libgtcore/cstr.h"
#include "libgtcore/ma.h"
#include "libgtext/mapping.h"

struct Mapping {
  Str *mapping_file;
  char *global;
  MappingType type;
  lua_State *L;
  bool is_table;
};

Mapping* mapping_new(Str *mapping_file, const char *global_name,
                     MappingType type, Error *e)
{
  Mapping *m;
  int had_err = 0;
  error_check(e);
  assert(mapping_file && global_name);
  /* alloc */
  m = ma_malloc(sizeof (Mapping));
  m->mapping_file = str_ref(mapping_file);
  m->global = cstr_dup(global_name);
  m->type = type;
  /* create new lua state (i.e., interpreter) */
  m->L = luaL_newstate();
  if (!m->L) {
    error_set(e, "out of memory (cannot create new lua state)");
    had_err = -1;
  }
  /* load the standard libs into the lua interpreter */
  if (!had_err)
    luaL_openlibs(m->L);
  /* try to load & run mapping file */
  if (!had_err) {
    if (luaL_loadfile(m->L, str_get(mapping_file)) ||
        lua_pcall(m->L, 0, 0, 0)) {
      error_set(e, "cannot run file: %s", lua_tostring(m->L, -1));
      had_err = -1;
    }
  }
  /* make sure a global variable with name <global_name> is defined */
  if (!had_err) {
    lua_getglobal(m->L, global_name);
    if (lua_isnil(m->L, -1)) {
      error_set(e, "'%s' is not defined in \"%s\"", global_name,
                str_get(mapping_file));
      had_err = -1;
    }
  }
  /* make sure it is either a table or a function */
  if (!had_err) {
    if (!(lua_istable(m->L, -1) || lua_isfunction(m->L, -1))) {
      error_set(e, "'%s' must be either a table or a function (defined "
                   "in \"%s\")", global_name, str_get(mapping_file));
      had_err = -1;
    }
  }
  /* remember if it is a table or a function */
  if (!had_err) {
    if (lua_istable(m->L, -1))
      m->is_table = true;
    else {
      m->is_table = false;
      lua_pop(m->L, 1);
    }
  }
  /* return */
  if (had_err) {
    mapping_delete(m);
    return NULL;
  }
  return m;
}

static int map_table(Mapping *m, Str **stroutput, long *integeroutput,
                     const char *input, Error *e)
{
  int had_err = 0;
  error_check(e);
  assert(m && input);
  assert((m->type == MAPPINGTYPE_STRING  && stroutput) ||
         (m->type == MAPPINGTYPE_INTEGER && integeroutput));
  lua_pushstring(m->L, input);
  lua_gettable(m->L, -2); /* get global[input] */
  /* make sure global[input] is defined */
  if (lua_isnil(m->L, -1)) {
    error_set(e, "%s[%s] is nil (defined in \"%s\")", m->global, input,
              str_get(m->mapping_file)); had_err = -1;
  }
  if (!had_err) {
    switch (m->type) {
      case MAPPINGTYPE_STRING:
        /* make sure global[input] is a string */
        if (!(lua_isstring(m->L, -1))) {
          error_set(e, "%s[%s] is not a string (defined in \"%s\")", m->global,
                    input, str_get(m->mapping_file)); had_err = -1;
        }
        if (!had_err)
          *stroutput = str_new_cstr(lua_tostring(m->L, -1));
        break;
      case MAPPINGTYPE_INTEGER:
        /* make sure global[input] is an integer */
        if (!(lua_isnumber(m->L, -1))) {
          error_set(e, "%s[%s] is not an integer (defined in \"%s\")",
                    m->global, input, str_get(m->mapping_file)); had_err = -1;
        }
        if (!had_err)
          *integeroutput = lua_tointeger(m->L, -1);
        break;
    }
  }
  lua_pop(m->L, 1); /* pop result */
  return had_err;
}

static int map_function(Mapping *m, Str **stroutput, long *integeroutput,
                        const char *input, Error *e)
{
  int had_err = 0;
  error_check(e);
  assert(m && input);
  assert((m->type == MAPPINGTYPE_STRING  && stroutput) ||
         (m->type == MAPPINGTYPE_INTEGER && integeroutput));
  lua_getglobal(m->L, m->global);
  lua_pushstring(m->L, input);
  /* call function */
  if (lua_pcall(m->L, 1, 1, 0)) {
    error_set(e, "running function '%s': %s", m->global,
              lua_tostring(m->L, -1));
    had_err = -1;
  }
  /* make sure the result is a string */
  if (!had_err) {
    switch (m->type) {
      case MAPPINGTYPE_STRING:
        if (!lua_isstring(m->L, -1)) {
          error_set(e, "function '%s' must return a string (defined in \"%s\")",
                    m->global, str_get(m->mapping_file));
          had_err = -1;
        }
        if (!had_err)
          *stroutput = str_new_cstr(lua_tostring(m->L, -1));
         break;
       case MAPPINGTYPE_INTEGER:
        if (!lua_isnumber(m->L, -1)) {
          error_set(e, "function '%s' must return an integer) (defined in "
                       "\"%s\")", m->global, str_get(m->mapping_file));
          had_err = -1;
        }
        if (!had_err)
          *integeroutput = lua_tointeger(m->L, -1);
        break;
    }
  }
  lua_pop(m->L, 1); /* pop result */
  return had_err;
}

static int map_generic(Mapping *m, Str **stroutput, long *integeroutput,
                       const char *input, Error *e)
{
  error_check(e);
  assert(m && input);
  assert((m->type == MAPPINGTYPE_STRING  && stroutput) ||
         (m->type == MAPPINGTYPE_INTEGER && integeroutput));
  if (m->is_table)
    return map_table(m, stroutput, integeroutput, input, e);
  return map_function(m, stroutput, integeroutput, input, e);
}

Str* mapping_map_string(Mapping *m, const char *input, Error *e)
{
  Str *output = NULL;
  error_check(e);
  map_generic(m, &output, NULL, input, e);
  return output;
}

int mapping_map_integer(Mapping *m, long *output, const char *input, Error *e)
{
  error_check(e);
  return map_generic(m, NULL, output, input, e);
}

void mapping_delete(Mapping *m)
{
  if (!m) return;
  str_delete(m->mapping_file);
  ma_free(m->global);
  if (m->L) lua_close(m->L);
  ma_free(m);
}
