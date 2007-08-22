/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GTLUA_H
#define GTLUA_H

#include "lua.h"

int  luaopen_gt(lua_State*); /* open all GenomeTools libraries in Lua */
void run_interactive_lua_interpreter(lua_State*);

#endif
