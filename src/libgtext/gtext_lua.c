/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "libgtext/genome_node_lua.h"
#include "libgtext/genome_stream_lua.h"
#include "libgtext/gtext_lua.h"

int luaopen_gtext(lua_State *L)
{
  assert(L);
  luaopen_genome_node(L);
  luaopen_genome_stream(L);
  return 1;
}
