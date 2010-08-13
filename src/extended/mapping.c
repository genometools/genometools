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

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "core/assert_api.h"
#include "core/cstr_api.h"
#include "core/ma.h"
#include "extended/mapping.h"

struct GtMapping {
  GtStr *mapping_file;
  char *global;
  GtMappingType type;
  lua_State *L;
  bool is_table;
};

GtMapping* gt_mapping_new(GtStr *mapping_file, const char *global_name,
                          GtMappingType type, GtError *err)
{
  GtMapping *m;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(mapping_file && global_name);
  /* alloc */
  m = gt_malloc(sizeof (GtMapping));
  m->mapping_file = gt_str_ref(mapping_file);
  m->global = gt_cstr_dup(global_name);
  m->type = type;
  /* create new lua state (i.e., interpreter) */
  m->L = luaL_newstate();
  if (!m->L) {
    gt_error_set(err, "out of memory (cannot create new Lua state)");
    had_err = -1;
  }
  /* load the standard libs into the Lua interpreter */
  if (!had_err) {
    luaL_openlibs(m->L);
    gt_assert(!lua_gettop(m->L));
  }
  /* try to load & run mapping file */
  if (!had_err) {
    if (luaL_loadfile(m->L, gt_str_get(mapping_file)) ||
        lua_pcall(m->L, 0, 0, 0)) {
      gt_error_set(err, "cannot run file: %s", lua_tostring(m->L, -1));
      had_err = -1;
      lua_pop(m->L, 1);
    }
  }
  /* make sure a global variable with name <global_name> is defined */
  if (!had_err) {
    lua_getglobal(m->L, global_name);
    if (lua_isnil(m->L, -1)) {
      gt_error_set(err, "'%s' is not defined in \"%s\"", global_name,
                gt_str_get(mapping_file));
      had_err = -1;
      lua_pop(m->L, 1);
    }
  }
  /* make sure it is either a table or a function */
  if (!had_err) {
    if (!(lua_istable(m->L, -1) || lua_isfunction(m->L, -1))) {
      gt_error_set(err, "'%s' must be either a table or a function (defined "
                     "in \"%s\")", global_name, gt_str_get(mapping_file));
      had_err = -1;
      lua_pop(m->L, 1);
    }
  }
  /* remember if it is a table or a function */
  if (!had_err) {
    if (lua_istable(m->L, -1))
      m->is_table = true;
    else
      m->is_table = false;
    lua_pop(m->L, 1);
  }
  /* return */
  gt_assert(!lua_gettop(m->L));
  if (had_err) {
    gt_mapping_delete(m);
    return NULL;
  }
  return m;
}

static int map_table(GtMapping *m, GtStr **stroutput, long *integeroutput,
                     const char *input, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(m && input);
  gt_assert((m->type == GT_MAPPINGTYPE_STRING  && stroutput) ||
            (m->type == GT_MAPPINGTYPE_INTEGER && integeroutput));
  gt_assert(!lua_gettop(m->L));
  lua_getglobal(m->L, m->global);
  lua_pushstring(m->L, input);
  lua_gettable(m->L, -2); /* get global[input] */
  /* make sure global[input] is defined */
  if (lua_isnil(m->L, -1)) {
    gt_error_set(err, "%s[%s] is nil (defined in \"%s\")", m->global, input,
                 gt_str_get(m->mapping_file));
    had_err = -1;
  }
  if (!had_err) {
    switch (m->type) {
      case GT_MAPPINGTYPE_STRING:
        /* make sure global[input] is a string */
        if (!(lua_isstring(m->L, -1))) {
          gt_error_set(err, "%s[%s] is not a string (defined in \"%s\")",
                        m->global, input, gt_str_get(m->mapping_file));
          had_err = -1;
        }
        if (!had_err)
          *stroutput = gt_str_new_cstr(lua_tostring(m->L, -1));
        break;
      case GT_MAPPINGTYPE_INTEGER:
        /* make sure global[input] is an integer */
        if (!(lua_isnumber(m->L, -1))) {
          gt_error_set(err, "%s[%s] is not an integer (defined in \"%s\")",
                       m->global, input, gt_str_get(m->mapping_file));
          had_err = -1;
        }
        if (!had_err)
          *integeroutput = lua_tointeger(m->L, -1);
        break;
    }
  }
  lua_pop(m->L, 1); /* pop result */
  lua_pop(m->L, 1); /* pop table */
  gt_assert(!lua_gettop(m->L));
  return had_err;
}

static int map_function(GtMapping *m, GtStr **stroutput, long *integeroutput,
                        const char *input, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(m && input);
  gt_assert((m->type == GT_MAPPINGTYPE_STRING  && stroutput) ||
            (m->type == GT_MAPPINGTYPE_INTEGER && integeroutput));
  gt_assert(!lua_gettop(m->L));
  lua_getglobal(m->L, m->global);
  lua_pushstring(m->L, input);
  /* call function */
  if (lua_pcall(m->L, 1, 1, 0)) {
    gt_error_set(err, "running function '%s': %s", m->global,
              lua_tostring(m->L, -1));
    had_err = -1;
  }
  /* make sure the result is a string */
  if (!had_err) {
    switch (m->type) {
      case GT_MAPPINGTYPE_STRING:
        if (!lua_isstring(m->L, -1)) {
          gt_error_set(err, "function '%s' must return a string (defined in "
                         "\"%s\")", m->global, gt_str_get(m->mapping_file));
          had_err = -1;
        }
        if (!had_err)
          *stroutput = gt_str_new_cstr(lua_tostring(m->L, -1));
         break;
       case GT_MAPPINGTYPE_INTEGER:
        if (!lua_isnumber(m->L, -1)) {
          gt_error_set(err, "function '%s' must return an integer) (defined in "
                       "\"%s\")", m->global, gt_str_get(m->mapping_file));
          had_err = -1;
        }
        if (!had_err)
          *integeroutput = lua_tointeger(m->L, -1);
        break;
    }
  }
  lua_pop(m->L, 1); /* pop result */
  gt_assert(!lua_gettop(m->L));
  return had_err;
}

static int map_generic(GtMapping *m, GtStr **stroutput, long *integeroutput,
                       const char *input, GtError *err)
{
  gt_error_check(err);
  gt_assert(m && input);
  gt_assert((m->type == GT_MAPPINGTYPE_STRING  && stroutput) ||
            (m->type == GT_MAPPINGTYPE_INTEGER && integeroutput));
  if (m->is_table)
    return map_table(m, stroutput, integeroutput, input, err);
  return map_function(m, stroutput, integeroutput, input, err);
}

GtStr* gt_mapping_map_string(GtMapping *m, const char *input, GtError *err)
{
  GtStr *output = NULL;
  gt_error_check(err);
  map_generic(m, &output, NULL, input, err);
  return output;
}

int gt_mapping_map_integer(GtMapping *m, long *output, const char *input,
                           GtError *err)
{
  gt_error_check(err);
  return map_generic(m, NULL, output, input, err);
}

void gt_mapping_delete(GtMapping *m)
{
  if (!m) return;
  gt_str_delete(m->mapping_file);
  gt_free(m->global);
  if (m->L) lua_close(m->L);
  gt_free(m);
}
