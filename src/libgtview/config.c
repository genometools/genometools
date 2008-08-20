  /*
  Copyright (c) 2007-2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c)      2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/warning.h"
#include "libgtext/luahelper.h"
#include "libgtext/luaserialize.h"
#include "libgtview/config.h"
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

struct Config
{
  lua_State *L;
  Str *filename;
  bool verbose;
};

static void config_lua_new_table(lua_State *L, const char *key)
{
  lua_pushstring(L, key);
  lua_newtable(L);
  lua_settable(L, -3);
}

static const luaL_Reg luasecurelibs[] = {
  /* Think very hard before adding additional Lua libraries to this list, it
     might compromise the security of web applications like GenomeViewer!
     Do not add the 'io', 'os', or 'debug' library under any circumstances! */
  {"", luaopen_base},
  {LUA_TABLIBNAME, luaopen_table},
  {LUA_STRLIBNAME, luaopen_string},
  {LUA_MATHLIBNAME, luaopen_math},
  {NULL, NULL}
};

static void luaL_opensecurelibs(lua_State *L)
{
  const luaL_Reg *lib = luasecurelibs;
  for (; lib->func; lib++) {
    lua_pushcfunction(L, lib->func);
    lua_pushstring(L, lib->name);
    lua_call(L, 1, 0);
  }
}

Config* config_new(bool verbose, Error *err)
{
  Config *cfg;
  error_check(err);
  cfg = ma_calloc(1, sizeof (Config));
  cfg->filename = NULL;
  cfg->verbose = verbose;
  cfg->L = luaL_newstate();
  if (!cfg->L) {
    error_set(err, "out of memory (cannot create new lua state)");
    ma_free(cfg);
    return NULL;
  }
  else
    luaL_opensecurelibs(cfg->L); /* do not replace with luaL_openlibs()! */
  return cfg;
}

Config* config_new_with_state(lua_State *L)
{
  Config *cfg;
  cfg = ma_calloc(1, sizeof (Config));
  cfg->L = L;
  return cfg;
}

int config_load_file(Config *cfg, Str *fn, Error *err)
{
  int had_err = 0;
  error_check(err);
  assert(cfg && cfg->L && fn);
  cfg->filename = str_ref(fn);
  if (config_get_verbose(cfg))
    fprintf(stderr, "Trying to load config file: %s...\n", str_get(fn));
  if (luaL_loadfile(cfg->L, str_get(fn)) ||
      lua_pcall(cfg->L, 0, 0, 0)) {
    error_set(err, "cannot run configuration file: %s",
              lua_tostring(cfg->L, -1));
    had_err = -1;
  }
  if (!had_err) {
    lua_getglobal(cfg->L, "config");
    if (lua_isnil(cfg->L, -1) || !lua_istable(cfg->L, -1)) {
      error_set(err, "'config' is not defined or not a table in \"%s\"",
                str_get(fn));
    }
    lua_pop(cfg->L, 1);
  }
  return had_err;
}

void config_reload(Config *cfg)
{
  int rval;
  assert(cfg && cfg->filename);
  rval = config_load_file(cfg, cfg->filename, NULL);
  assert(!rval); /* should not happen, config file was loaded before */
}

/* Searches for <section> inside the config table, creating it if it does not
   exist and finally pushing it on the Lua stack (at the top).
   Returns the total number of items pushed on the stack by this function. */
static int config_find_section_for_setting(Config* cfg, const char *section)
{
  int depth = 0;
  assert(cfg && section);
  lua_getglobal(cfg->L, "config");
  if (lua_isnil(cfg->L, -1)) {
    lua_pop(cfg->L, 1);
    lua_newtable(cfg->L);
    lua_setglobal(cfg->L, "config");
    lua_getglobal(cfg->L, "config");
  }
  depth++;
  lua_getfield(cfg->L, -1, section);
  if (lua_isnil(cfg->L, -1)) {
    lua_pop(cfg->L, 1);
    config_lua_new_table(cfg->L, section);
    lua_getfield(cfg->L, -1, section);
  }
  depth++;
  return depth;
}

/* Searches for <section> inside the config table, returning -1 if it is not
   found. Otherwise the number of items pushed onto the stack is returned. */
static int config_find_section_for_getting(const Config *cfg,
                                           const char *section)
{
  int depth = 0;
  assert(cfg && section);
  lua_getglobal(cfg->L, "config");
  if (lua_isnil(cfg->L, -1)) {
    if (cfg->verbose) warning("'config' is not defined");
    lua_pop(cfg->L, 1);
    return -1;
  } else depth++;
  lua_getfield(cfg->L, -1, section);
  if (lua_isnil(cfg->L, -1) || !lua_istable(cfg->L, -1)) {
    if (cfg->verbose) warning("section '%s' is not defined", section);
    lua_pop(cfg->L, 1);
    return -1;
  } else depth++;
  return depth;
}

bool config_get_color(const Config *cfg, const char *section,
                      const char *key, Color *color)
{
  int i = 0;
  assert(cfg && section && key && color);
  /* set default colors */
  color->red=0.5; color->green = 0.5; color->blue=0.5;
  /* get section */
  i = config_find_section_for_getting(cfg, section);
  /* could not get section, return default */
  if (i < 0) {
    lua_pop(cfg->L, i);
    return false;
  }
  /* lookup color entry for given feature */
  lua_getfield(cfg->L, -1, key);
  if (lua_isnil(cfg->L, -1) || !lua_istable(cfg->L, -1)) {
    if (cfg->verbose) warning("no colors are defined for type '%s', "
                               "will use defaults.",
                               key);
    lua_pop(cfg->L, 3);
    return false;
  } else i++;
  /* update color struct */
  lua_getfield(cfg->L, -1, "red");
  if (lua_isnil(cfg->L, -1) || !lua_isnumber(cfg->L, -1)) {
    if (cfg->verbose) warning("%s  value for type '%s' is undefined or"
                               " not numeric, using default","red", key);
  }
  else
    color->red = lua_tonumber(cfg->L,-1);
  lua_pop(cfg->L, 1);
  lua_getfield(cfg->L, -1, "green");
  if (lua_isnil(cfg->L, -1) || !lua_isnumber(cfg->L, -1)) {
    if (cfg->verbose) warning("%s  value for type '%s' is undefined or"
                               " not numeric, using default","green", key);
  }
  else
    color->green = lua_tonumber(cfg->L,-1);
  lua_pop(cfg->L, 1);
  lua_getfield(cfg->L, -1, "blue");
  if (lua_isnil(cfg->L, -1) || !lua_isnumber(cfg->L, -1)) {
    if (cfg->verbose) warning("%s  value for type '%s' is undefined or"
                               " not numeric, using default","blue", key);
  }
  else
    color->blue = lua_tonumber(cfg->L,-1);
  lua_pop(cfg->L, 1);
  /* reset stack to original state for subsequent calls */
  lua_pop(cfg->L, i);
  return true;
}

void config_set_color(Config *cfg, const char *section, const char *key,
                      Color *color)
{
  int i = 0;
  assert(cfg && section && key && color);
  i = config_find_section_for_setting(cfg, section);
  lua_getfield(cfg->L, -1, key);
  i++;
  if (lua_isnil(cfg->L, -1)) {
    lua_pop(cfg->L, 1);
    config_lua_new_table(cfg->L, key);
    lua_getfield(cfg->L, -1, key);
  }
  lua_pushstring(cfg->L, "red");
  lua_pushnumber(cfg->L, color->red);
  lua_settable(cfg->L, -3);
  lua_pushstring(cfg->L, "green");
  lua_pushnumber(cfg->L, color->green);
  lua_settable(cfg->L, -3);
  lua_pushstring(cfg->L, "blue");
  lua_pushnumber(cfg->L, color->blue);
  lua_settable(cfg->L, -3);
  lua_pop(cfg->L, i);
}

bool config_get_str(const Config *cfg, const char *section,
                     const char *key, Str *text)
{
  int i = 0;
  assert(cfg && key && section);
  /* get section */
  i = config_find_section_for_getting(cfg, section);
  /* could not get section, return default */
  if (i < 0) {
    lua_pop(cfg->L, i);
    return false;
  }
  /* lookup entry for given key */
  lua_getfield(cfg->L, -1, key);
  if (lua_isnil(cfg->L, -1) || !lua_isstring(cfg->L, -1)) {
    if (cfg->verbose) warning("no value is defined for key '%s'",
                               key);
    lua_pop(cfg->L, i+1);
    return false;
  } else i++;
  /* retrieve string */
  str_set(text, lua_tostring(cfg->L, -1));
  /* reset stack to original state for subsequent calls */
  lua_pop(cfg->L, i);
  return true;
}

void config_set_str(Config *cfg, const char *section, const char *key,
                     Str *str)
{
  int i = 0;
  assert(cfg && section && key && str);
  i = config_find_section_for_setting(cfg, section);
  lua_pushstring(cfg->L, key);
  lua_pushstring(cfg->L, str_get(str));
  lua_settable(cfg->L, -3);
  lua_pop(cfg->L, i);
}

bool config_get_num(const Config *cfg, const char *section, const char *key,
                    double *val)
{
  int i = 0;
  assert(cfg && key && section && val);
  /* get section */
  i = config_find_section_for_getting(cfg, section);
  /* could not get section, return default */
  if (i < 0) {
    lua_pop(cfg->L, i);
    return false;
  }
  /* lookup entry for given key */
  lua_getfield(cfg->L, -1, key);
  if (lua_isnil(cfg->L, -1) || !lua_isnumber(cfg->L, -1)) {
    if (cfg->verbose) warning("no or non-numeric value found for key '%s'",
                              key);
    lua_pop(cfg->L, i+1);
    return false;
  } else i++;
  /* retrieve value */
  *val = lua_tonumber(cfg->L, -1);
  /* reset stack to original state for subsequent calls */
  lua_pop(cfg->L, i);
  return true;
}

void config_set_num(Config *cfg, const char *section, const char *key,
                    double number)
{
  int i = 0;
  assert(cfg && section && key);
  i = config_find_section_for_setting(cfg, section);
  lua_pushstring(cfg->L, key);
  lua_pushnumber(cfg->L, number);
  lua_settable(cfg->L, -3);
  lua_pop(cfg->L, i);
}

bool config_get_bool(const Config *cfg, const char *section, const char *key,
                     bool *val)
{
  int i = 0;
  assert(cfg && key && section);
  /* get section */
  i = config_find_section_for_getting(cfg, section);
  /* could not get section, return default */
  if (i < 0) {
    lua_pop(cfg->L, i);
    return false;
  }
  /* lookup entry for given key */
  lua_getfield(cfg->L, -1, key);
  if (lua_isnil(cfg->L, -1) || !lua_isboolean(cfg->L, -1)) {
    if (cfg->verbose) warning("no or non-boolean value found for key '%s'",
                              key);
    lua_pop(cfg->L, i+1);
    return false;
  } else i++;
  /* retrieve value */
  *val = lua_toboolean(cfg->L, -1);
  /* reset stack to original state for subsequent calls */
  lua_pop(cfg->L, i);
  return true;
}

void config_set_bool(Config *cfg, const char *section, const char *key,
                     bool flag)
{
  int i = 0;
  assert(cfg && section && key);
  i = config_find_section_for_setting(cfg, section);
  lua_pushstring(cfg->L, key);
  lua_pushboolean(cfg->L, flag);
  lua_settable(cfg->L, -3);
  lua_pop(cfg->L, i);
}

void config_unset(Config *cfg, const char *section, const char *key)
{
  assert(cfg && section && key);
  lua_getglobal(cfg->L, "config");
  if (!lua_isnil(cfg->L, -1)) {
    assert(lua_istable(cfg->L, -1));
    lua_getfield(cfg->L, -1, section);
    if (!lua_isnil(cfg->L, -1)) {
      assert(lua_istable(cfg->L, -1));
      lua_pushstring(cfg->L, key);
      lua_pushnil(cfg->L);
      lua_settable(cfg->L, -3);
    }
  }
}

bool config_get_verbose(const Config *cfg)
{
  assert(cfg);
  return cfg->verbose;
}

int config_to_str(const Config *cfg, Str *outstr)
{
  int had_err = 0;
  assert(cfg && outstr);
  lua_getglobal(cfg->L, "config");
  str_append_cstr(outstr, "config = {\n");
  had_err = lua_table_to_str(cfg->L, outstr, -1);
  str_append_cstr(outstr, "}");
  lua_pop(cfg->L, 1);
  return had_err;
}

int config_load_str(Config *cfg, Str *instr)
{
  int had_err = 0;
  assert(cfg && instr);
  had_err = luaL_loadbuffer(cfg->L, str_get(instr), str_length(instr), "str") ||
              lua_pcall(cfg->L, 0, 0, 0);
  return had_err;
}

Config* config_clone(const Config *cfg, Error *err)
{
  int had_err = 0;
  Str *cfg_buffer = str_new();
  Config *new_cfg = NULL;
  assert(cfg);
  if (!(new_cfg = config_new(config_get_verbose(cfg), err)))
    had_err = -1;
  if (!had_err)
    had_err = config_to_str(cfg, cfg_buffer);
  if (!had_err)
    had_err = config_load_str(new_cfg, cfg_buffer);
  if (had_err)
      error_set(err, "An error occurred trying to clone Config object %p", cfg);
  str_delete(cfg_buffer);
  return new_cfg;
}

int config_unit_test(Error *err)
{
  int had_err = 0;
  Config *cfg = NULL, *new_cfg = NULL;
  bool val;
  Str *luafile = str_new_cstr("config.lua"),
      *test1   = str_new_cstr("mRNA"),
      *str     = str_new(),
      *cfg_buffer = str_new();
  Color col1, col2, col, defcol, tmpcol;
  double num;
  error_check(err);

  /* example colors */
  col1.red=.1;col1.green=.2;col1.blue=.3;
  col2.red=.4;col2.green=.5;col2.blue=.6;
  col.red=1.0;col.green=1.0;col.blue=1.0;
  defcol.red=.5;defcol.green=.5;defcol.blue=.5;

  /* instantiate new config object */
  if (!(cfg = config_new(false, err)))
    had_err = -1;

  /* at the beginning, all values are defaults, since nothing is defined */
  config_get_color(cfg, "exon", "fill", &tmpcol);
  ensure(had_err, color_equals(tmpcol,defcol));
  config_get_color(cfg, "cds", "fill", &tmpcol);
  ensure(had_err, color_equals(tmpcol,defcol));
  config_get_color(cfg, "foo", "fill", &tmpcol);
  ensure(had_err, color_equals(tmpcol,defcol));
  if (!config_get_num(cfg, "format", "margins", &num))
    num = 10.0;
  ensure(had_err, num == 10.0);
  if (!config_get_bool(cfg, "exon", "collapse_to_parent", &val))
    val = false;
  ensure(had_err, !val);

  /* change some values... */
  config_set_color(cfg, "exon", "fill", &col1);
  config_set_num(cfg, "format", "margins", 11.0);
  config_set_num(cfg, "format", "foo", 2.0);

  /* is it saved correctly? */
  config_get_color(cfg, "exon", "fill", &tmpcol);
  ensure(had_err, !color_equals(tmpcol,defcol));
  ensure(had_err, color_equals(tmpcol,col1));
  if (!config_get_num(cfg, "format", "margins", &num))
    num = 10.0;
  ensure(had_err, num == 11.0);
  if (!config_get_num(cfg, "format", "foo", &num))
    num = 2.0;
  ensure(had_err, num == 2.0);

  /* create a new color definition */
  config_set_color(cfg, "foo", "fill", &col2);
  config_set_str(cfg, "bar", "baz", test1);

  /* is it saved correctly? */
  config_get_color(cfg, "foo", "fill", &tmpcol);
  ensure(had_err, !color_equals(tmpcol,defcol));
  ensure(had_err, color_equals(tmpcol,col2));
  if (!config_get_str(cfg, "bar", "baz", str))
    str_set(str, "");
  ensure(had_err, (strcmp(str_get(str),"")!=0));
  ensure(had_err, (str_cmp(str,test1)==0));
  if (!config_get_str(cfg, "bar", "test", str))
    str_set(str, "");
  ensure(had_err, (strcmp(str_get(str),"")==0));

  /* clone a Config object */
  new_cfg = config_clone(cfg, err);
  error_check(err);
  if (!error_is_set(err))
  {
    /* check again */
    config_get_color(new_cfg, "foo", "fill", &tmpcol);
    ensure(had_err, !color_equals(tmpcol,defcol));
    ensure(had_err, color_equals(tmpcol,col2));
    if (!config_get_str(new_cfg, "bar", "baz", str))
      str_set(str, "");
    ensure(had_err, (strcmp(str_get(str),"")!=0));
    ensure(had_err, (str_cmp(str,test1)==0));
    if (!config_get_str(new_cfg, "bar", "test", str))
      str_set(str, "");
    ensure(had_err, (strcmp(str_get(str),"")==0));
  }
  /* mem cleanup */
  str_delete(luafile);
  str_delete(test1);
  str_delete(str);
  str_delete(cfg_buffer);
  config_delete(cfg);
  config_delete(new_cfg);

  return had_err;
}

void config_delete_without_state(Config *cfg)
{
  if (!cfg) return;
  str_delete(cfg->filename);
  ma_free(cfg);
}

void config_delete(Config *cfg)
{
  if (!cfg) return;
  if (cfg->L) lua_close(cfg->L);
  config_delete_without_state(cfg);
}
