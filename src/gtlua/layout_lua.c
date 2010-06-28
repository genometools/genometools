/*
  Copyright (c) 2008-2010 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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

#ifndef WITHOUT_CAIRO

#include "lauxlib.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/layout.h"
#include "annotationsketch/luastyle.h"
#include "core/error.h"
#include "extended/luahelper.h"
#include "gtlua/canvas_lua.h"
#include "gtlua/diagram_lua.h"
#include "gtlua/layout_lua.h"

static int layout_lua_new(lua_State *L)
{
  GtLayout **layout;
  GtDiagram **diagram;
  unsigned int width;
  GtStyle *style;
  GtError *err;
  diagram = check_diagram(L, 1);
  width = luaL_checkint(L, 2);
  /* create layout */
  style = gt_lua_get_style_from_registry(L);
  layout = lua_newuserdata(L, sizeof (GtLayout*));
  gt_assert(layout);
  err = gt_error_new();
  *layout = gt_layout_new(*diagram, width, style, err);
  if (gt_error_is_set(err))
    return gt_lua_error(L, err);
  gt_error_delete(err);
  luaL_getmetatable(L, LAYOUT_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int layout_lua_sketch(lua_State *L)
{
  GtLayout **layout;
  GtCanvas **canvas;
  GtError *err;
  int had_err = 0;
  layout = check_layout(L, 1);
  canvas = check_canvas(L, 2);
  err = gt_error_new();
  had_err = gt_layout_sketch(*layout, *canvas, err);
  if (had_err < 0)
    return gt_lua_error(L, err);
  gt_error_delete(err);
  return 0;
}

static int layout_lua_get_height(lua_State *L)
{
  GtLayout **layout;
  GtError *err;
  int had_err = 0;
  unsigned long height;
  layout = check_layout(L, 1);
  err = gt_error_new();
  had_err = gt_layout_get_height(*layout, &height, err);
  lua_pushnumber(L, height);
  if (had_err < 0)
    return gt_lua_error(L, err);
  gt_error_delete(err);
  return 1;
}

static int layout_lua_delete(lua_State *L)
{
  GtLayout **layout;
  layout = check_layout(L, 1);
  gt_layout_delete(*layout);
  return 0;
}

static const struct luaL_Reg layout_lib_f [] = {
  { "layout_new", layout_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg layout_lib_m [] = {
  { "sketch", layout_lua_sketch },
  { "get_height", layout_lua_get_height },
  { NULL, NULL }
};

int gt_lua_open_layout(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_newmetatable(L, LAYOUT_METATABLE);
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, layout_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, layout_lib_m);
  lua_pop(L, 1);
  luaL_register(L, "gt", layout_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}

#endif
