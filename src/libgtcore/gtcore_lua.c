/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "libgtcore/bittab_lua.h"
#include "libgtcore/gtcore_lua.h"

int luaopen_gtcore(lua_State *L)
{
  assert(L);
  luaopen_bittab(L);
  return 1;
}
