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

#ifndef HELPER_H
#define HELPER_H

#include "lua.h"
#include "libgtcore/env.h"

#ifdef LIBGTVIEW
#include "libgtview/config.h"

void put_config_in_registry(lua_State*, Config*);
Config* get_config_from_registry(lua_State*);
#endif

void put_env_in_registry(lua_State*, Env*);
Env* get_env_from_registry(lua_State*);
void set_arg_in_lua_interpreter(lua_State*, const char *argv_0,
                                const char **argv);
void run_interactive_lua_interpreter(lua_State*);
/* Propagate the error given in <env> (which must be set) to <L>. The error in
   <env> is unset. */
int  luagt_error(lua_State *L, Env *env);

#endif
