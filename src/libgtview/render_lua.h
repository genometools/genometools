/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef RENDER_LUA_H
#define RENDER_LUA_H

#include "lua.h"

/* exports the Render class to Lua:

   render = gt.render_new()
            render:to_png(diagram, filename, [width])
*/
int luaopen_render(lua_State*);

#endif
