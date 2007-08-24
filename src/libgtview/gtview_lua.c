/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "libgtview/feature_index_lua.h"
#include "libgtview/gtview_lua.h"

int luaopen_gtext(lua_State *L)
{
  assert(L);
  luaopen_feature_index(L);
  return 1;
}
