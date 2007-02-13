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
  Str *mapping_filename;
  lua_State *L;
  bool is_table;
};

RegionMapping* regionmapping_new(Str *mapping_filename)
{
  RegionMapping *rm;
  assert(mapping_filename);
  rm = xmalloc(sizeof(RegionMapping));
  rm->mapping_filename = str_ref(mapping_filename);
  rm->L = luaL_newstate(); /* create new lua state (i.e., interpreter) */
  if (!rm->L)
    error("out of memory (cannot create new lua state)");
  luaL_openlibs(rm->L); /* load the standard libs into the lua interpreter */
  /* try to load & run mapping file */
  if (luaL_loadfile(rm->L, str_get(mapping_filename)) ||
      lua_pcall(rm->L, 0, 0, 0)) {
    error("cannot run mapping file: %s", lua_tostring(rm->L, -1));
  }
  lua_getglobal(rm->L, "mapping");
  /* make sure a global 'mapping' variable is defined */
  if (lua_isnil(rm->L, -1))
    error("'mapping' is not defined in \"%s\"", str_get(mapping_filename));
  /* make sure it is either a table or a function */
  if (!(lua_istable(rm->L, -1) || lua_isfunction(rm->L, -1))) {
    error("'mapping' must be either a table or a function (defined in \"%s\")",
          str_get(mapping_filename));
  }
  if (lua_istable(rm->L, -1))
    rm->is_table = true;
  else
    rm->is_table = false;
  return rm;
}

Str* regionmapping_map(RegionMapping *rm, const char *sequence_region)
{
  Str *result;
  assert(rm && sequence_region);
  if (rm->is_table) {
    lua_pushstring(rm->L, sequence_region);
    lua_gettable(rm->L, -2); /* get mapping[sequence_region] */
    /* make sure mapping[sequence_region] is defined */
    if (lua_isnil(rm->L, -1)) {
      error("mapping[%s] is nil (defined in \"%s\")", sequence_region,
            str_get(rm->mapping_filename));
    }
    /* make sure mapping[sequence_region] is a string */
    if (!(lua_isstring(rm->L, -1))) {
      error("mapping[%s] is not a string (defined in \"%s\")", sequence_region,
            str_get(rm->mapping_filename));
    }
    result = str_new_cstr(lua_tostring(rm->L, -1));
    lua_pop(rm->L, 1);
  }
  else {
    result = NULL; /* XXX */
  }
  return result;
}

void regionmapping_free(RegionMapping *rm)
{
  if (!rm) return;
  str_free(rm->mapping_filename);
  lua_close(rm->L);
  free(rm);
}
