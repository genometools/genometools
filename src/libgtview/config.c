  /*
  Copyright (c) 2007      Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
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

Config* config_new(bool verbose, Error *err)
{
  Config *cfg;
  error_check(err);
  cfg = ma_calloc(1, sizeof (Config));
  cfg->filename = NULL;
  cfg->verbose = verbose;
  cfg->L = luaL_newstate();
  if (!cfg->L)
  {
    error_set(err, "out of memory (cannot create new lua state)");
  } else luaL_openlibs(cfg->L);
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

Color config_get_color(const Config *cfg, const char *key)
{
  Color color;
  int i = 0;
  assert(cfg && key);
  /* set default colors */
  color.red=0.8; color.green = 0.8; color.blue=0.8;
  /* get section */
  i = config_find_section_for_getting(cfg, "colors");
  /* could not get section, return default */
  if (i < 0) {
    lua_pop(cfg->L, i);
    return color;
  }
  /* lookup color entry for given feature */
  lua_getfield(cfg->L, -1, key);
  if (lua_isnil(cfg->L, -1) || !lua_istable(cfg->L, -1)) {
    if (cfg->verbose) warning("no colors are defined for type '%s', "
                               "will use defaults.",
                               key);
    lua_pop(cfg->L, 3);
    return color;
  } else i++;
  /* update color struct */
  lua_getfield(cfg->L, -1, "red");
  if (lua_isnil(cfg->L, -1) || !lua_isnumber(cfg->L, -1)) {
    if (cfg->verbose) warning("%s  value for type '%s' is undefined or"
                               " not numeric, using default","red", key);
  }
  else
    color.red = lua_tonumber(cfg->L,-1);
  lua_pop(cfg->L, 1);
  lua_getfield(cfg->L, -1, "green");
  if (lua_isnil(cfg->L, -1) || !lua_isnumber(cfg->L, -1)) {
    if (cfg->verbose) warning("%s  value for type '%s' is undefined or"
                               " not numeric, using default","green", key);
  }
  else
    color.green = lua_tonumber(cfg->L,-1);
  lua_pop(cfg->L, 1);
  lua_getfield(cfg->L, -1, "blue");
  if (lua_isnil(cfg->L, -1) || !lua_isnumber(cfg->L, -1)) {
    if (cfg->verbose) warning("%s  value for type '%s' is undefined or"
                               " not numeric, using default","blue", key);
  }
  else
    color.blue = lua_tonumber(cfg->L,-1);
  lua_pop(cfg->L, 1);
  /* reset stack to original state for subsequent calls */
  lua_pop(cfg->L, i);
  return color;
}

void config_get_colorptr(const Config *cfg, Color *color, const char *key)
{
  assert(cfg && color && key);
  *color = config_get_color(cfg, key);
}

void config_set_color(Config *cfg, const char *key, Color *color)
{
  int i = 0;
  assert(cfg && key);
  i = config_find_section_for_setting(cfg, "colors");
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

const char* config_get_cstr(const Config *cfg, const char *section,
                            const char *key, const char *deflt)
{
  const char *str = deflt;
  int i = 0;
  assert(cfg && key && section);
  /* get section */
  i = config_find_section_for_getting(cfg, section);
  /* could not get section, return default */
  if (i < 0) {
    lua_pop(cfg->L, i);
    return deflt;
  }
  /* lookup entry for given key */
  lua_getfield(cfg->L, -1, key);
  if (lua_isnil(cfg->L, -1) || !lua_isstring(cfg->L, -1)) {
    if (cfg->verbose) warning("no value is defined for key '%s'",
                               key);
    lua_pop(cfg->L, i+1);
    return deflt;
  } else i++;
  /* retrieve string */
  str = lua_tostring(cfg->L, -1);
  /* reset stack to original state for subsequent calls */
  lua_pop(cfg->L, i);
  return str;
}

void config_set_cstr(Config *cfg, const char *section, const char *key,
                     const char *str)
{
  int i = 0;
  assert(cfg && section && key);
  i = config_find_section_for_setting(cfg, section);
  lua_pushstring(cfg->L, key);
  lua_pushstring(cfg->L, str);
  lua_settable(cfg->L, -3);
  lua_pop(cfg->L, i);
}

double config_get_num(const Config *cfg, const char *section, const char *key,
                      double deflt)
{
  double num = deflt;
  int i = 0;
  assert(cfg && key && section);
  /* get section */
  i = config_find_section_for_getting(cfg, section);
  /* could not get section, return default */
  if (i < 0) {
    lua_pop(cfg->L, i);
    return deflt;
  }
  /* lookup entry for given key */
  lua_getfield(cfg->L, -1, key);
  if (lua_isnil(cfg->L, -1) || !lua_isnumber(cfg->L, -1)) {
    if (cfg->verbose) warning("no or non-numeric value found for key '%s'",
                              key);
    lua_pop(cfg->L, i+1);
    return deflt;
  } else i++;
  /* retrieve value */
  num = lua_tonumber(cfg->L, -1);
  /* reset stack to original state for subsequent calls */
  lua_pop(cfg->L, i);
  return num;
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

StrArray* config_get_cstr_list(const Config *cfg, const char *section,
                               const char *key)
{
  StrArray *list = NULL;
  int i, had_err = 0;
  assert(cfg && section && key);
  i = config_find_section_for_getting(cfg, section);
  if (i < 0)
    had_err = -1;
  if (!had_err) {
    /* lookup entry for given key */
    lua_getfield(cfg->L, -1, key);
    if (!lua_istable(cfg->L, -1))
      had_err = -1;
  }
  /* build list */
  if (!had_err) {
    lua_Integer i = 0;
    list = strarray_new();
    for (;;) {
      lua_pushinteger(cfg->L, ++i);
      lua_gettable(cfg->L, -2);
      if (lua_isnil(cfg->L, -1))
        break;
      if (lua_isstring(cfg->L, -1))
        strarray_add_cstr(list, lua_tostring(cfg->L, -1));
      lua_pop(cfg->L, 1);
    }
  }
  /* return result */
  if (had_err) {
    strarray_delete(list);
    return NULL;
  }
  return list;
}

void config_set_cstr_list(Config *cfg,  const char *section, const char *key,
                          StrArray *list)
{
  assert(cfg && section && key && list);
  config_find_section_for_setting(cfg, section);
  assert(lua_istable(cfg->L, -1));
  lua_pushstring(cfg->L, key);
  lua_push_strarray_as_table(cfg->L, list);
  lua_settable(cfg->L, -3);
}

bool config_cstr_in_list(const Config *cfg, const char *section,
                         const char *key, const char *checkstr)
{
  int i = 0, had_err = 0;
  bool ret = false;
  assert(cfg && key && section && checkstr);
  /* get section */
  i = config_find_section_for_getting(cfg, section);
  if (i < 0) {
    lua_pop(cfg->L, i);
    return false;
  }
  lua_getfield(cfg->L, -1, key);
  i++;
  if (lua_isnil(cfg->L, -1) || !lua_istable(cfg->L, -1)) {
    lua_pop(cfg->L, 1);
    i--;
    had_err = -1;
  }
  if (!had_err)
  {
    /* table is at the top position in the stack */
    lua_pushnil(cfg->L);  /* first key */
    while (lua_next(cfg->L, -2) != 0) {
      /* uses 'key' (at index -2) and 'value' (at index -1) */
      if (lua_isnil(cfg->L, -1) || !lua_isstring(cfg->L, -1))
      {
        if (cfg->verbose) warning("non-string value in section %s, key %s",
                                   section, key);
        lua_pop(cfg->L, 1);
      } else {
        if (strcmp(checkstr, lua_tostring(cfg->L, -1)) == 0)
        {
          /* found checkstr, report success */
          ret = true;
        }
        /* removes 'value'; keeps 'key' for next iteration */
        lua_pop(cfg->L, 1);
      }
    }
  }
  /* remove temporary stack items */
  lua_pop(cfg->L, i);
  return ret;
}

bool config_get_verbose(const Config *cfg)
{
  assert(cfg);
  return cfg->verbose;
}

DominateStatus config_dominates(Config* cfg, GenomeFeatureType gft1,
                                GenomeFeatureType gft2)
{
  char *fts1, *fts2;

  assert(cfg);

  if (gft1 == gft2)
    return DOMINATES_EQUAL;
  fts1 = (char*) genome_feature_type_get_cstr(gft1);
  fts2 = (char*) genome_feature_type_get_cstr(gft2);
  if (fts1 == NULL || fts2 == NULL)
    return DOMINATES_UNKNOWN_TYPE;
  else {
    if (config_cstr_in_list(cfg, "dominate",fts1, fts2))
      return DOMINATES_FIRST;
    else if (config_cstr_in_list(cfg, "dominate",fts2, fts1))
      return DOMINATES_SECOND;
    else
      return DOMINATES_NOT_SPECIFIED;
  }
}

int config_unit_test(Error *err)
{
  int had_err = 0;
  Config *cfg;
  const char* test1 = "mRNA";
  const char* str = NULL;
  Str *luafile = str_new_cstr("config.lua");
  Color col1, col2, col, defcol, tmpcol;
  double num;
  error_check(err);

  /* example colors */
  col1.red=.1;col1.green=.2;col1.blue=.3;
  col2.red=.4;col2.green=.5;col2.blue=.6;
  col.red=1.0;col.green=1.0;col.blue=1.0;
  defcol.red=.8;defcol.green=.8;defcol.blue=.8;

  /* instantiate new config object */
  if (!(cfg = config_new(false, err)))
    had_err = -1;

  /* at the beginning, all values are defaults, since nothing is defined */
  tmpcol = config_get_color(cfg, "exon");
  ensure(had_err, color_equals(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "cds");
  ensure(had_err, color_equals(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "foo");
  ensure(had_err, color_equals(tmpcol,defcol));
  num = config_get_num(cfg,"format", "margins", 10.0);
  ensure(had_err, num == 10.0);
  str = config_get_cstr(cfg, "collapse", "exon", "");
  ensure(had_err, (strcmp(str,"")==0));

  /* change some values... */
  config_set_color(cfg, "exon", &col);
  config_set_num(cfg,"format", "margins", 11.0);
  config_set_num(cfg,"format", "foo", 2.0);

  /* is it saved correctly? */
  tmpcol = config_get_color(cfg, "exon");
  ensure(had_err, !color_equals(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "exon");
  ensure(had_err, color_equals(tmpcol,col));
  num = config_get_num(cfg,"format", "margins", 10.0);
  ensure(had_err, num == 11.0);
  num = config_get_num(cfg,"format", "foo", 10.0);
  ensure(had_err, num == 2.0);

  /* create a new color definition */
  config_set_color(cfg, "foo", &col);
  config_set_cstr(cfg, "bar", "baz", test1);

  /* is it saved correctly? */
  tmpcol = config_get_color(cfg, "foo");
  ensure(had_err, !color_equals(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "foo");
  ensure(had_err, color_equals(tmpcol,col));
  str = config_get_cstr(cfg, "bar", "baz", "");
  ensure(had_err, (strcmp(str,"")!=0));
  ensure(had_err, (strcmp(str,test1)==0));
  str = config_get_cstr(cfg, "bar", "test", "");
  ensure(had_err, (strcmp(str,"")==0));

  /* mem cleanup */
  str_delete(luafile);
  config_delete(cfg);

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
