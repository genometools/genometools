/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/symbol.h"
#include "core/unused_api.h"
#include "extended/script_filter.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/gt_lua.h"
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

struct GtScriptFilter
{
  lua_State *L;
  GtStr *filename;
  unsigned long reference_count;
};

static const luaL_Reg script_filter_luasecurelibs[] = {
  {"", luaopen_base},
  {LUA_TABLIBNAME, luaopen_table},
  {LUA_STRLIBNAME, luaopen_string},
  {LUA_MATHLIBNAME, luaopen_math},
  {"gt", gt_lua_open_lib}, /* open the GenomeTools library for callbacks */
  {NULL, NULL}
};

static void script_filter_luaL_opencustomlibs(lua_State *L,
                                              const luaL_Reg *lib)
{
  for (; lib->func; lib++) {
    lua_pushcfunction(L, lib->func);
    lua_pushstring(L, lib->name);
    lua_call(L, 1, 0);
  }
}

GtScriptFilter *gt_script_filter_new(const char *file, GtError *err)
{
  GtScriptFilter *script_filter;
  gt_assert(file);
  script_filter = gt_malloc(sizeof (GtScriptFilter));
  script_filter->filename = gt_str_new_cstr(file);
  script_filter->L = luaL_newstate();
  script_filter->reference_count = 0;
  if (!script_filter->L) {
    gt_error_set(err, "out of memory (cannot create new Lua state)");
    gt_free(script_filter);
    return NULL;
  }
  script_filter_luaL_opencustomlibs(script_filter->L,
                                     script_filter_luasecurelibs);
  if (luaL_loadfile(script_filter->L, file) ||
                               lua_pcall(script_filter->L, 0, 0, 0)) {
    gt_error_set(err, "cannot run file: %s",
                 lua_tostring(script_filter->L, -1));
    lua_pop(script_filter->L, 1);
    lua_close(script_filter->L);
    gt_str_delete(script_filter->filename);
    gt_free(script_filter);
    return NULL;
  }
  return script_filter;
}

GtScriptFilter *gt_script_filter_new_from_string(const char *script_string,
                                                 GtError *err)
{
  GtScriptFilter *script_filter;
  gt_assert(script_string);
  script_filter = gt_malloc(sizeof (GtScriptFilter));
  script_filter->filename = NULL;
  script_filter->L = luaL_newstate();
  script_filter->reference_count = 0;
  if (!script_filter->L) {
    gt_error_set(err, "out of memory (cannot create new Lua state)");
    gt_free(script_filter);
    return NULL;
  }
  script_filter_luaL_opencustomlibs(script_filter->L,
                                     script_filter_luasecurelibs);
  if (luaL_loadstring(script_filter->L, script_string) ||
                               lua_pcall(script_filter->L, 0, 0, 0)) {
    gt_error_set(err, "cannot run file: %s",
                 lua_tostring(script_filter->L, -1));
    lua_pop(script_filter->L, 1);
    lua_close(script_filter->L);
    gt_free(script_filter);
    return NULL;
  }
  return script_filter;
}

/* TODO: caching */
static const char *gt_script_filter_get_string(GtScriptFilter *script_filter,
                                              const char *name, GtError *err)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(script_filter && name);
  gt_error_check(err);

#ifndef NDEBUG
  stack_size = lua_gettop(script_filter->L);
#endif

  lua_getglobal(script_filter->L, name);
  if (lua_isnil(script_filter->L, -1)) {
    lua_pop(script_filter->L, 1);
    return gt_symbol("undefined");
  }
  /* execute callback if function is given */
  if (lua_isfunction(script_filter->L, -1))
  {
    int num_of_args = 0;
    if (lua_pcall(script_filter->L, num_of_args, 1, 0) != 0)
    {
      gt_error_set(err, "%s", lua_tostring(script_filter->L, -1));
      lua_pop(script_filter->L, 1);
      gt_assert(lua_gettop(script_filter->L) == stack_size);
      return NULL;
    }
  }

  if (lua_isnil(script_filter->L, -1) || !lua_isstring(script_filter->L, -1)) {
    lua_pop(script_filter->L, 1);
    gt_assert(lua_gettop(script_filter->L) == stack_size);
    gt_error_set(err, "script filter '%s': '%s' must return a string",
                 gt_str_get(script_filter->filename), name);
    return NULL;
  }

  /* retrieve string */
  return lua_tostring(script_filter->L, -1);
}

const char* gt_script_filter_get_name(GtScriptFilter *script_filter,
                                      GtError *err)
{
  return gt_script_filter_get_string(script_filter, "name", err);
}

const char* gt_script_filter_get_description(GtScriptFilter *script_filter,
                                             GtError *err)
{
  return gt_script_filter_get_string(script_filter, "description", err);
}

const char*
gt_script_filter_get_short_description(GtScriptFilter *script_filter,
                                       GtError *err)
{
  return gt_script_filter_get_string(script_filter, "short_descr", err);
}

const char* gt_script_filter_get_author(GtScriptFilter *script_filter,
                                        GtError *err)
{
  return gt_script_filter_get_string(script_filter, "author", err);
}

const char* gt_script_filter_get_email(GtScriptFilter *script_filter,
                                       GtError *err)
{
  return gt_script_filter_get_string(script_filter, "email", err);
}

const char* gt_script_filter_get_version(GtScriptFilter *script_filter,
                                         GtError *err)
{
  return gt_script_filter_get_string(script_filter, "version", err);
}

bool gt_script_filter_validate(GtScriptFilter *script_filter, GtError *err)
{
  const char *result;

#ifndef NDEBUG
  GT_UNUSED int stack_size;
#endif
  gt_assert(script_filter);
  gt_error_check(err);

#ifndef NDEBUG
  stack_size = lua_gettop(script_filter->L);
#endif

  result = gt_script_filter_get_name(script_filter, err);
  if (result == gt_symbol("undefined")) {
    gt_error_set(err, "metadata 'name' not found");
    return false;
  }
  result = gt_script_filter_get_description(script_filter, err);
  if (result == gt_symbol("undefined")) {
    gt_error_set(err, "metadata 'description' not found");
    return false;
  }
  result = gt_script_filter_get_short_description(script_filter, err);
  if (result == gt_symbol("undefined")) {
    gt_error_set(err, "metadata 'short_descr' not found");
    return false;
  }
  result = gt_script_filter_get_author(script_filter, err);
  if (result == gt_symbol("undefined")) {
    gt_error_set(err, "metadata 'author' not found");
    return false;
  }
  result = gt_script_filter_get_email(script_filter, err);
  if (result == gt_symbol("undefined")) {
    gt_error_set(err, "metadata 'email' not found");
    return false;
  }
  result = gt_script_filter_get_version(script_filter, err);
  if (result == gt_symbol("undefined")) {
    gt_error_set(err, "metadata 'version' not found");
    return false;
  }

  lua_getglobal(script_filter->L, "filter");
  if (lua_isnil(script_filter->L, -1)) {
    gt_error_set(err, "function 'filter' is not defined");
    lua_pop(script_filter->L, 1);
    return false;
  }
  return true;
}

int gt_script_filter_run(GtScriptFilter *sf, GtFeatureNode *gf,
                         bool *select_node, GtError *err)
{
  int had_err = 0;
#ifndef NDEBUG
  int stack_size;
#endif
  GtGenomeNode *gn_lua;

#ifndef NDEBUG
  stack_size = lua_gettop(sf->L);
#endif

  if (!had_err) {
    lua_getglobal(sf->L, "filter");
    if (lua_isnil(sf->L, -1)) {
      gt_error_set(err, "function 'filter' is not defined");
      had_err = -1;
      lua_pop(sf->L, 1);
    }
  }

  if (!had_err) {
    gn_lua = gt_genome_node_ref((GtGenomeNode*) gf);
    gt_lua_genome_node_push(sf->L, gn_lua);

    if (lua_pcall(sf->L, 1, 1, 0) != 0) {
      gt_error_set(err, "error running function 'filter': %s",
                   lua_tostring(sf->L, -1));
      lua_pop(sf->L, 1);
      had_err = -1;
    }
  }

  if (!had_err && !lua_isboolean(sf->L, -1)) {
    gt_error_set(err, "function 'filter' must return boolean");
    lua_pop(sf->L, 1);
    had_err = -1;
  }

  if (!had_err) {
    *select_node = lua_toboolean(sf->L, -1);
    lua_pop(sf->L, 1);
  }
  gt_assert(lua_gettop(sf->L) == stack_size);

  return had_err;
}

GtScriptFilter* gt_script_filter_ref(GtScriptFilter *script_filter)
{
  if (!script_filter) return NULL;
  script_filter->reference_count++;
  return script_filter;
}

void gt_script_filter_delete(GtScriptFilter *script_filter)
{
  if (!script_filter) return;
  if (script_filter->reference_count) {
    script_filter->reference_count--;
    return;
  }
  gt_str_delete(script_filter->filename);
  lua_close(script_filter->L);
  gt_free(script_filter);
}
