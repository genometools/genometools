/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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

#ifndef WITHOUT_CAIRO

#include "lauxlib.h"
#include "annotationsketch/canvas_cairo_file.h"
#include "annotationsketch/luastyle.h"
#include "core/error.h"
#include "extended/luahelper.h"
#include "gtlua/canvas_lua.h"
#include "gtlua/image_info_lua.h"

static int canvas_cairo_file_lua_new_generic(lua_State *L, GtGraphicsOutType t)
{
  GtCanvas **canvas;
  GtImageInfo **ii;
  unsigned int width,
               height;
  GtError *err;
  GtStyle *style;
  width = luaL_checkint(L, 1);
  height = luaL_checkint(L, 2);
  /* create canvas */
  style = gt_lua_get_style_from_registry(L);
  canvas = lua_newuserdata(L, sizeof (GtCanvas*));
  gt_assert(canvas);
  /* if a imageinfo object is passed, it must be correct type */
  if (lua_isnil(L, 3)) {
    err = gt_error_new();
    *canvas = gt_canvas_cairo_file_new(style, t, width, height, NULL, err);
  } else {
    ii = check_imageinfo(L, 3);
    err = gt_error_new();
    *canvas = gt_canvas_cairo_file_new(style, t, width, height, *ii, err);
  }
  if (gt_error_is_set(err))
    return gt_lua_error(L, err);
  gt_error_delete(err);
  luaL_getmetatable(L, CANVAS_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int canvas_cairo_file_lua_new_pdf(lua_State *L)
{
  return canvas_cairo_file_lua_new_generic(L, GT_GRAPHICS_PDF);
}

static int canvas_cairo_file_lua_new_png(lua_State *L)
{
  return canvas_cairo_file_lua_new_generic(L, GT_GRAPHICS_PNG);
}

static int canvas_cairo_file_lua_new_svg(lua_State *L)
{
  return canvas_cairo_file_lua_new_generic(L, GT_GRAPHICS_SVG);
}

static int canvas_cairo_file_lua_new_ps(lua_State *L)
{
  return canvas_cairo_file_lua_new_generic(L, GT_GRAPHICS_PS);
}

static int canvas_cairo_file_lua_to_file(lua_State *L)
{
  GtCanvas **canvas;
  GtCanvasCairoFile *ccf = NULL;
  GtError *err;
  const char *fn;
  int had_err = 0;
  canvas = check_canvas(L, 1);
  ccf = canvas_cairo_file_try_cast(*canvas);
  luaL_argcheck(L, ccf, 1, "must be a CanvasCairoFile object");
  fn = luaL_checkstring(L, 2);
  gt_assert(canvas);
  err = gt_error_new();
  had_err = gt_canvas_cairo_file_to_file(ccf, fn, err);
  if (had_err)
    return gt_lua_error(L, err);
  gt_error_delete(err);
  return 0;
}

static int canvas_lua_delete(lua_State *L)
{
  GtCanvas **canvas;
  canvas = check_canvas(L, 1);
  gt_canvas_delete(*canvas);
  return 0;
}

static const struct luaL_Reg canvas_lib_f [] = {
  { "canvas_cairo_file_new_png", canvas_cairo_file_lua_new_png },
  { "canvas_cairo_file_new_pdf", canvas_cairo_file_lua_new_pdf },
  { "canvas_cairo_file_new_ps", canvas_cairo_file_lua_new_ps },
  { "canvas_cairo_file_new_svg", canvas_cairo_file_lua_new_svg },
  { NULL, NULL }
};

static const struct luaL_Reg canvas_lib_m [] = {
  { "to_file", canvas_cairo_file_lua_to_file },
  { NULL, NULL }
};

int gt_lua_open_canvas(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_newmetatable(L, CANVAS_METATABLE);
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, canvas_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, canvas_lib_m);
  lua_pop(L, 1);
  luaL_register(L, "gt", canvas_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}

#endif
