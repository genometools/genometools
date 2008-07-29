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

#ifdef LIBGTVIEW

#include "lauxlib.h"
#include "libgtcore/error.h"
#include "libgtext/luahelper.h"
#include "libgtlua/canvas_lua.h"
#include "libgtview/canvas.h"
#include "libgtview/luaconfig.h"

static int canvas_lua_new_png(lua_State *L)
{
  Canvas **canvas;
  /* ImageInfo **ii; */
  unsigned int width;
  Config *config;
  width = luaL_checkint(L, 1);
  /* create canvas */
  config = lua_get_config_from_registry(L);
  canvas = lua_newuserdata(L, sizeof (Canvas*));
  assert(canvas);
  /* TODO: ImageInfo */
  *canvas = canvas_new(config, GRAPHICS_PNG, width, NULL);
  luaL_getmetatable(L, CANVAS_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int canvas_lua_to_file(lua_State *L)
{
  Canvas **canvas;
  Error *err;
  const char *fn;
  int had_err = 0;
  err = error_new();
  canvas = check_canvas(L, 1);
  fn = luaL_checkstring(L, 2);
  assert(canvas);
  had_err = canvas_to_file(*canvas, fn, err);
  if (had_err)
    return lua_gt_error(L, err);
  error_delete(err);
  return 0;
}

static int canvas_lua_delete(lua_State *L)
{
  Canvas **canvas;
  canvas = check_canvas(L, 1);
  canvas_delete(*canvas);
  return 0;
}

static const struct luaL_Reg canvas_lib_f [] = {
  { "canvas_new_png", canvas_lua_new_png },
  { NULL, NULL }
};

static const struct luaL_Reg canvas_lib_m [] = {
  { "to_file", canvas_lua_to_file },
  { NULL, NULL }
};

int luaopen_canvas(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, CANVAS_METATABLE);
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, canvas_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, canvas_lib_m);
  luaL_register(L, "gt", canvas_lib_f);
  return 1;
}

#endif
