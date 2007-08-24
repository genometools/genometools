/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GTLUA_H
#define GTLUA_H

#include "lua.h"
#include "libgtcore/env.h"

#ifdef LIBGTVIEW
#include "libgtview/config.h"

void put_config_in_registry(lua_State*, Config*);
Config* get_config_from_registry(lua_State*);
#endif

void put_env_in_registry(lua_State*, Env*);
Env* get_env_from_registry(lua_State*);
int  luaopen_gt(lua_State*); /* open all GenomeTools libraries in Lua */
void set_arg_in_lua_interpreter(lua_State*, const char *argv_0,
                                const char **argv);
void run_interactive_lua_interpreter(lua_State*);
/* Propagate the error given in <env> (which must be set) to <L>. The error in
   <env> is unset. */
int  luagt_error(lua_State *L, Env *env);

#endif
