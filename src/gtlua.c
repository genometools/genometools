/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "gtlua.h"
#include "lauxlib.h"
#include "libgtcore/gtcorelua.h"

int luaopen_gt(lua_State *L)
{
  assert(L);
  luaopen_gtcore(L); /* open core library */
  return 1;
}

void run_interactive_lua_interpreter(lua_State *L)
{
  char buf[BUFSIZ];
  int error;
  assert(L);
  while (fgets(buf, sizeof buf, stdin)) {
    error = luaL_loadbuffer(L, buf, strlen(buf), "line") ||
            lua_pcall(L, 0, 0, 0);
    if (error) {
      fprintf(stderr, "%s", lua_tostring(L, -1));
      lua_pop(L, 1); /* pop error message */
    }
  }
}
