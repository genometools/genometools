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

#ifdef LIBGTVIEW

#include "lauxlib.h"
#include "libgtext/luahelper.h"
#include "libgtlua/diagram_lua.h"
#include "libgtlua/render_lua.h"
#include "libgtview/luaconfig.h"
#include "libgtview/render.h"

#define RENDER_METATABLE  "GenomeTools.render"
#define check_render(L) \
        (Render**) luaL_checkudata(L, 1, RENDER_METATABLE)

static int render_lua_new(lua_State *L)
{
  Render **render;
  Config *config;
  render = lua_newuserdata(L, sizeof (Render*));
  assert(render);
  config = lua_get_config_from_registry(L);
  *render = render_new(config);
  luaL_getmetatable(L, RENDER_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int render_lua_to_png(lua_State *L)
{
  Render **render;
  Diagram **diagram;
  const char *filename;
  unsigned int width;
  Error *err;
  render = check_render(L);
  diagram = check_diagram(L, 2);
  filename = luaL_checkstring(L, 3);
  if (lua_gettop(L) >= 4)
    width = luaL_checkint(L, 4);
  else
    width = DEFAULT_RENDER_WIDTH;
  err = error_new();
  if (render_to_png(*render, *diagram, filename, width, err))
    return lua_gt_error(L, err);
  error_delete(err);
  return 0;
}

static int render_lua_delete(lua_State *L)
{
  Render **render = check_render(L);
  render_delete(*render);
  return 0;
}

static const struct luaL_Reg render_lib_f [] = {
  { "render_new", render_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg render_lib_m [] = {
  { "to_png", render_lua_to_png },
  { NULL, NULL }
};

int luaopen_render(lua_State *L)
{
  assert(L);
  luaL_newmetatable(L, RENDER_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, render_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, render_lib_m);
  luaL_register(L, "gt", render_lib_f);
  return 1;
}

#endif
