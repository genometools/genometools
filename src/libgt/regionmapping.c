/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "error.h"
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "regionmapping.h"
#include "xansi.h"

struct RegionMapping {
  lua_State *L;
};

RegionMapping* regionmapping_new(const char *mapping_filename)
{
  RegionMapping *rm;
  assert(mapping_filename);
  rm = xmalloc(sizeof(RegionMapping));
  rm->L = luaL_newstate(); /* create new lua state (i.e., interpreter) */
  if (!rm->L)
    error("out of memory (cannot create new lua state)");
  luaL_openlibs(rm->L); /* load the standard libs into the lua interpreter */
  return rm;
}

const char* regionmapping_map(RegionMapping *rm, const char *sequence_region)
{
  assert(rm && sequence_region);
  assert(0); /* XXX */
  return NULL; /* XXX */
}

void regionmapping_free(RegionMapping *rm)
{
  if (!rm) return;
  lua_close(rm->L);
  free(rm);
}
