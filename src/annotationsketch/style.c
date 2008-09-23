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
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "annotationsketch/color.h"
#include "annotationsketch/style.h"
#include "core/cstr.h"
#include "core/ensure.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/luahelper.h"
#include "extended/luaserialize.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/gt_lua.h"

struct GtStyle
{
  lua_State *L;
  char *filename;
};

static void gt_style_lua_new_table(lua_State *L, const char *key)
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
  {"gt", gt_lua_open_lib}, /* open the GenomeTools library for callbacks */
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

GtStyle* gt_style_new(GtError *err)
{
  GtStyle *sty;
  gt_error_check(err);
  sty = gt_calloc(1, sizeof (GtStyle));
  sty->filename = NULL;
  sty->L = luaL_newstate();
  if (!sty->L) {
    gt_error_set(err, "out of memory (cannot create new Lua state)");
    gt_free(sty);
    return NULL;
  }
  else
    luaL_opensecurelibs(sty->L); /* do not replace with luaL_openlibs()! */
  return sty;
}

GtStyle* gt_style_new_with_state(lua_State *L)
{
  GtStyle *sty;
  sty = gt_calloc(1, sizeof (GtStyle));
  sty->L = L;
  return sty;
}

int gt_style_load_file(GtStyle *sty, const char *filename, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  assert(sty && sty->L && filename);
  sty->filename = gt_cstr_dup(filename);
  gt_log_log("Trying to load style file: %s...", filename);
  if (luaL_loadfile(sty->L, filename) || lua_pcall(sty->L, 0, 0, 0)) {
    gt_error_set(err, "cannot run style file: %s",
              lua_tostring(sty->L, -1));
    had_err = -1;
  }
  if (!had_err) {
    lua_getglobal(sty->L, "style");
    if (lua_isnil(sty->L, -1) || !lua_istable(sty->L, -1)) {
      gt_error_set(err, "'style' is not defined or not a table in \"%s\"",
                filename);
    }
    lua_pop(sty->L, 1);
  }
  return had_err;
}

void gt_style_reload(GtStyle *sty)
{
  int rval;
  assert(sty && sty->filename);
  rval = gt_style_load_file(sty, sty->filename, NULL);
  assert(!rval); /* should not happen, file was loaded before */
}

/* Searches for <section> inside the style table, creating it if it does not
   exist and finally pushing it on the Lua stack (at the top).
   Returns the total number of items pushed on the stack by this function. */
static int gt_style_find_section_for_setting(GtStyle* sty, const char *section)
{
  int depth = 0;
  assert(sty && section);
  lua_getglobal(sty->L, "style");
  if (lua_isnil(sty->L, -1)) {
    lua_pop(sty->L, 1);
    lua_newtable(sty->L);
    lua_setglobal(sty->L, "style");
    lua_getglobal(sty->L, "style");
  }
  depth++;
  lua_getfield(sty->L, -1, section);
  if (lua_isnil(sty->L, -1)) {
    lua_pop(sty->L, 1);
    gt_style_lua_new_table(sty->L, section);
    lua_getfield(sty->L, -1, section);
  }
  depth++;
  return depth;
}

/* Searches for <section> inside the style table, returning -1 if it is not
   found. Otherwise the number of items pushed onto the stack is returned. */
static int gt_style_find_section_for_getting(const GtStyle *sty,
                                          const char *section)
{
  int depth = 0;
  assert(sty && section);
  lua_getglobal(sty->L, "style");
  if (lua_isnil(sty->L, -1)) {
    gt_warning("'style' is not defined");
    lua_pop(sty->L, 1);
    return -1;
  } else depth++;
  lua_getfield(sty->L, -1, section);
  if (lua_isnil(sty->L, -1) || !lua_istable(sty->L, -1)) {
    gt_warning("section '%s' is not defined", section);
    lua_pop(sty->L, 2);
    return -1;
  } else depth++;
  return depth;
}

bool gt_style_get_color(const GtStyle *sty, const char *section,
                     const char *key, GtColor *color, GtFeatureNode *gn)
{
  int i = 0;
  assert(sty && section && key && color);
  /* set default colors */
  color->red=0.5; color->green = 0.5; color->blue=0.5;
  /* get section */
  i = gt_style_find_section_for_getting(sty, section);
  /* could not get section, return default */
  if (i < 0) {
    return false;
  }
  /* lookup color entry for given feature */
  lua_getfield(sty->L, -1, key);

  /* execute callback if function is given */
  if (lua_isfunction(sty->L, -1) && gn)
  {
    GtGenomeNode *gn_lua = gt_genome_node_rec_ref((GtGenomeNode*) gn);
    gt_lua_genome_node_push(sty->L, gn_lua);
    if (lua_pcall(sty->L, 1, 1, 0) != 0)
    {
      lua_pop(sty->L, 3);
      return false;
    }
  }

  if (lua_isnil(sty->L, -1) || !lua_istable(sty->L, -1)) {
    gt_warning("no colors are defined for type '%s', will use defaults.", key);
    lua_pop(sty->L, 3);
    return false;
  } else i++;
  /* update color struct */
  lua_getfield(sty->L, -1, "red");
  if (lua_isnil(sty->L, -1) || !lua_isnumber(sty->L, -1)) {
    gt_warning("%s  value for type '%s' is undefined or not numeric, using "
               "default","red", key);
  }
  else
    color->red = lua_tonumber(sty->L,-1);
  lua_pop(sty->L, 1);
  lua_getfield(sty->L, -1, "green");
  if (lua_isnil(sty->L, -1) || !lua_isnumber(sty->L, -1)) {
    gt_warning("%s  value for type '%s' is undefined or not numeric, using "
               "default","green", key);
  }
  else
    color->green = lua_tonumber(sty->L,-1);
  lua_pop(sty->L, 1);
  lua_getfield(sty->L, -1, "blue");
  if (lua_isnil(sty->L, -1) || !lua_isnumber(sty->L, -1)) {
    gt_warning("%s  value for type '%s' is undefined or not numeric, using "
               "default","blue", key);
  }
  else
    color->blue = lua_tonumber(sty->L,-1);
  lua_pop(sty->L, 1);
  /* reset stack to original state for subsequent calls */
  lua_pop(sty->L, i);
  return true;
}

void gt_style_set_color(GtStyle *sty, const char *section, const char *key,
                        const GtColor *color)
{
  int i = 0;
  assert(sty && section && key && color);
  i = gt_style_find_section_for_setting(sty, section);
  lua_getfield(sty->L, -1, key);
  i++;
  if (lua_isnil(sty->L, -1)) {
    lua_pop(sty->L, 1);
    gt_style_lua_new_table(sty->L, key);
    lua_getfield(sty->L, -1, key);
  }
  lua_pushstring(sty->L, "red");
  lua_pushnumber(sty->L, color->red);
  lua_settable(sty->L, -3);
  lua_pushstring(sty->L, "green");
  lua_pushnumber(sty->L, color->green);
  lua_settable(sty->L, -3);
  lua_pushstring(sty->L, "blue");
  lua_pushnumber(sty->L, color->blue);
  lua_settable(sty->L, -3);
  lua_pop(sty->L, i);
}

bool gt_style_get_str(const GtStyle *sty, const char *section,
                     const char *key, GtStr *text, GtFeatureNode *gn)
{
  int i = 0;
  assert(sty && key && section);
  /* get section */
  i = gt_style_find_section_for_getting(sty, section);
  /* could not get section, return default */
  if (i < 0) {
    return false;
  }
  /* lookup entry for given key */
  lua_getfield(sty->L, -1, key);

  /* execute callback if function is given */
  if (lua_isfunction(sty->L, -1) && gn)
  {
    GtGenomeNode *gn_lua = gt_genome_node_rec_ref((GtGenomeNode*) gn);
    gt_lua_genome_node_push(sty->L, gn_lua);
    if (lua_pcall(sty->L, 1, 1, 0) != 0)
    {
      lua_pop(sty->L, 3);
      return false;
    }
  }

  if (lua_isnil(sty->L, -1) || !lua_isstring(sty->L, -1)) {
    gt_warning("no value is defined for key '%s'", key);
    lua_pop(sty->L, i+1);
    return false;
  } else i++;
  /* retrieve string */
  gt_str_set(text, lua_tostring(sty->L, -1));
  /* reset stack to original state for subsequent calls */
  lua_pop(sty->L, i);
  return true;
}

void gt_style_set_str(GtStyle *sty, const char *section, const char *key,
                     GtStr *str)
{
  int i = 0;
  assert(sty && section && key && str);
  i = gt_style_find_section_for_setting(sty, section);
  lua_pushstring(sty->L, key);
  lua_pushstring(sty->L, gt_str_get(str));
  lua_settable(sty->L, -3);
  lua_pop(sty->L, i);
}

bool gt_style_get_num(const GtStyle *sty, const char *section, const char *key,
                    double *val, GT_UNUSED GtFeatureNode *gn)
{
  int i = 0;
  assert(sty && key && section && val);
  /* get section */
  i = gt_style_find_section_for_getting(sty, section);
  /* could not get section, return default */
  if (i < 0) {
    return false;
  }
  /* lookup entry for given key */
  lua_getfield(sty->L, -1, key);

  /* execute callback if function is given */
  if (lua_isfunction(sty->L, -1) && gn)
  {
    GtGenomeNode *gn_lua = gt_genome_node_rec_ref((GtGenomeNode*) gn);
    gt_lua_genome_node_push(sty->L, gn_lua);
    if (lua_pcall(sty->L, 1, 1, 0) != 0)
    {
      lua_pop(sty->L, 3);
      return false;
    }
  }

  if (lua_isnil(sty->L, -1) || !lua_isnumber(sty->L, -1)) {
    gt_warning("no or non-numeric value found for key '%s'", key);
    lua_pop(sty->L, i+1);
    return false;
  } else i++;
  /* retrieve value */
  *val = lua_tonumber(sty->L, -1);
  /* reset stack to original state for subsequent calls */
  lua_pop(sty->L, i);
  return true;
}

void gt_style_set_num(GtStyle *sty, const char *section, const char *key,
                    double number)
{
  int i = 0;
  assert(sty && section && key);
  i = gt_style_find_section_for_setting(sty, section);
  lua_pushstring(sty->L, key);
  lua_pushnumber(sty->L, number);
  lua_settable(sty->L, -3);
  lua_pop(sty->L, i);
}

bool gt_style_get_bool(const GtStyle *sty, const char *section,
                       const char *key, bool *val,
                       GT_UNUSED GtFeatureNode *gn)
{
  int i = 0;
  assert(sty && key && section);
  /* get section */
  i = gt_style_find_section_for_getting(sty, section);
  /* could not get section, return default */
  if (i < 0) {
    return false;
  }
  /* lookup entry for given key */
  lua_getfield(sty->L, -1, key);
  if (lua_isnil(sty->L, -1) || !lua_isboolean(sty->L, -1)) {
    gt_warning("no or non-boolean value found for key '%s'", key);
    lua_pop(sty->L, i+1);
    return false;
  } else i++;
  /* retrieve value */
  *val = lua_toboolean(sty->L, -1);
  /* reset stack to original state for subsequent calls */
  lua_pop(sty->L, i);
  return true;
}

void gt_style_set_bool(GtStyle *sty, const char *section, const char *key,
                     bool flag)
{
  int i = 0;
  assert(sty && section && key);
  i = gt_style_find_section_for_setting(sty, section);
  lua_pushstring(sty->L, key);
  lua_pushboolean(sty->L, flag);
  lua_settable(sty->L, -3);
  lua_pop(sty->L, i);
}

void gt_style_unset(GtStyle *sty, const char *section, const char *key)
{
  assert(sty && section && key);
  lua_getglobal(sty->L, "style");
  if (!lua_isnil(sty->L, -1)) {
    assert(lua_istable(sty->L, -1));
    lua_getfield(sty->L, -1, section);
    if (!lua_isnil(sty->L, -1)) {
      assert(lua_istable(sty->L, -1));
      lua_pushstring(sty->L, key);
      lua_pushnil(sty->L);
      lua_settable(sty->L, -3);
    }
  }
}

int gt_style_to_str(const GtStyle *sty, GtStr *outstr, GtError *err)
{
  int had_err;
  gt_error_check(err);
  assert(sty && outstr);
  lua_getglobal(sty->L, "style");
  gt_str_append_cstr(outstr, "style = {\n");
  if (lua_istable(sty->L, -1))
    had_err = gt_lua_table_to_str(sty->L, outstr, -1, err);
  else {
    gt_error_set(err, "'style' must be a table");
    had_err = -1;
  }
  gt_str_append_cstr(outstr, "}");
  lua_pop(sty->L, 1);
  return had_err;
}

int gt_style_load_str(GtStyle *sty, GtStr *instr, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  assert(sty && instr);
  if (luaL_loadbuffer(sty->L, gt_str_get(instr), gt_str_length(instr), "str") ||
      lua_pcall(sty->L, 0, 0, 0)) {
    gt_error_set(err, "cannot run style buffer: %s",
              lua_tostring(sty->L, -1));
    had_err = -1;
  }
  return had_err;
}

GtStyle* gt_style_clone(const GtStyle *sty, GtError *err)
{
  int had_err = 0;
  GtStr *sty_buffer = gt_str_new();
  GtStyle *new_sty;
  assert(sty);
  if (!(new_sty = gt_style_new(err)))
    had_err = -1;
  if (!had_err)
    had_err = gt_style_to_str(sty, sty_buffer, err);
  if (!had_err)
    had_err = gt_style_load_str(new_sty, sty_buffer, err);
  gt_str_delete(sty_buffer);
  return new_sty;
}

int gt_style_unit_test(GtError *err)
{
  int had_err = 0;
  GtStyle *sty = NULL, *new_sty = NULL;
  bool val;
  GtStr *test1   = gt_str_new_cstr("mRNA"),
      *str     = gt_str_new(),
      *sty_buffer = gt_str_new();
  GtColor col1, col2, col, defcol, tmpcol;
  double num;
  gt_error_check(err);

  /* example colors */
  col1.red=.1;col1.green=.2;col1.blue=.3;
  col2.red=.4;col2.green=.5;col2.blue=.6;
  col.red=1.0;col.green=1.0;col.blue=1.0;
  defcol.red=.5;defcol.green=.5;defcol.blue=.5;

  /* instantiate new style object */
  if (!(sty = gt_style_new(err)))
    had_err = -1;

  /* at the beginning, all values are defaults, since nothing is defined */
  gt_style_get_color(sty, "exon", "fill", &tmpcol, NULL);
  ensure(had_err, gt_color_equals(&tmpcol, &defcol));
  gt_style_get_color(sty, "cds", "fill", &tmpcol, NULL);
  ensure(had_err, gt_color_equals(&tmpcol, &defcol));
  gt_style_get_color(sty, "foo", "fill", &tmpcol, NULL);
  ensure(had_err, gt_color_equals(&tmpcol, &defcol));
  if (!gt_style_get_num(sty, "format", "margins", &num, NULL))
    num = 10.0;
  ensure(had_err, num == 10.0);
  if (!gt_style_get_bool(sty, "exon", "collapse_to_parent", &val, NULL))
    val = false;
  ensure(had_err, !val);

  /* change some values... */
  gt_style_set_color(sty, "exon", "fill", &col1);
  gt_style_set_num(sty, "format", "margins", 11.0);
  gt_style_set_num(sty, "format", "foo", 2.0);

  /* is it saved correctly? */
  gt_style_get_color(sty, "exon", "fill", &tmpcol, NULL);
  ensure(had_err, !gt_color_equals(&tmpcol, &defcol));
  ensure(had_err, gt_color_equals(&tmpcol, &col1));
  if (!gt_style_get_num(sty, "format", "margins", &num, NULL))
    num = 10.0;
  ensure(had_err, num == 11.0);
  if (!gt_style_get_num(sty, "format", "foo", &num, NULL))
    num = 2.0;
  ensure(had_err, num == 2.0);

  /* create a new color definition */
  gt_style_set_color(sty, "foo", "fill", &col2);
  gt_style_set_str(sty, "bar", "baz", test1);

  /* is it saved correctly? */
  gt_style_get_color(sty, "foo", "fill", &tmpcol, NULL);
  ensure(had_err, !gt_color_equals(&tmpcol, &defcol));
  ensure(had_err, gt_color_equals(&tmpcol, &col2));
  if (!gt_style_get_str(sty, "bar", "baz", str, NULL))
    gt_str_set(str, "");
  ensure(had_err, (strcmp(gt_str_get(str),"")!=0));
  ensure(had_err, (gt_str_cmp(str,test1)==0));
  if (!gt_style_get_str(sty, "bar", "test", str, NULL))
    gt_str_set(str, "");
  ensure(had_err, (strcmp(gt_str_get(str),"")==0));

  /* clone a GtStyle object */
  new_sty = gt_style_clone(sty, err);
  gt_error_check(err);
  if (!gt_error_is_set(err))
  {
    /* check again */
    gt_style_get_color(new_sty, "foo", "fill", &tmpcol, NULL);
    ensure(had_err, !gt_color_equals(&tmpcol, &defcol));
    ensure(had_err, gt_color_equals(&tmpcol, &col2));
    if (!gt_style_get_str(new_sty, "bar", "baz", str, NULL))
      gt_str_set(str, "");
    ensure(had_err, (strcmp(gt_str_get(str),"")!=0));
    ensure(had_err, (gt_str_cmp(str,test1)==0));
    if (!gt_style_get_str(new_sty, "bar", "test", str, NULL))
      gt_str_set(str, "");
    ensure(had_err, (strcmp(gt_str_get(str),"")==0));
  }
  /* mem cleanup */
  gt_str_delete(test1);
  gt_str_delete(str);
  gt_str_delete(sty_buffer);
  gt_style_delete(sty);
  gt_style_delete(new_sty);

  return had_err;
}

void gt_style_delete_without_state(GtStyle *sty)
{
  if (!sty) return;
  gt_free(sty->filename);
  gt_free(sty);
}

void gt_style_delete(GtStyle *sty)
{
  if (!sty) return;
  if (sty->L) lua_close(sty->L);
  gt_style_delete_without_state(sty);
}
