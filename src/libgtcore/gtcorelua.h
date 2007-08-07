/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GTCORELUA_H
#define GTCORELUA_H

#include "lua.h"
#include "libgtcore/env.h"

void put_env_in_registry(lua_State*, Env*);
Env* get_env_from_registry(lua_State*);
int  luaopen_gtcore(lua_State*); /* open core library in Lua */

#endif
