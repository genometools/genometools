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

RegionMapping* regionmapping_new(Str *mapping_filename, Error *err)
{
  RegionMapping *rm;
  int has_err = 0;
  error_check(err);
  assert(mapping_filename);
  /* alloc */
  rm = xcalloc(1, sizeof (RegionMapping));
  rm->mapping_filename = str_ref(mapping_filename);
  /* create new lua state (i.e., interpreter) */
  rm->L = luaL_newstate();
  if (!rm->L) {
    /* XXX: -> assert? */
    error_set(err, "out of memory (cannot create new lua state)");
    has_err = -1;
  }
  /* load the standard libs into the lua interpreter */
  if (!has_err)
    luaL_openlibs(rm->L);
  /* try to load & run mapping file */
  if (!has_err) {
    if (luaL_loadfile(rm->L, str_get(mapping_filename)) ||
        lua_pcall(rm->L, 0, 0, 0)) {
      error_set(err, "cannot run mapping file: %s", lua_tostring(rm->L, -1));
      has_err = -1;
    }
  }
  /* make sure a global 'mapping' variable is defined */
  if (!has_err) {
    lua_getglobal(rm->L, "mapping");
    if (lua_isnil(rm->L, -1)) {
      error_set(err, "'mapping' is not defined in \"%s\"",
                str_get(mapping_filename));
      has_err = -1;
    }
  }
  /* make sure it is either a table or a function */
  if (!has_err) {
    if (!(lua_istable(rm->L, -1) || lua_isfunction(rm->L, -1))) {
      error_set(err, "'mapping' must be either a table or a function (defined "
                "in \"%s\")", str_get(mapping_filename));
      has_err = -1;
    }
  }
  /* remember if it is a table or a function */
  if (!has_err) {
    if (lua_istable(rm->L, -1))
      rm->is_table = true;
    else {
      rm->is_table = false;
      lua_pop(rm->L, 1);
    }
  }
  /* return */
  if (has_err) {
    regionmapping_free(rm);
    return NULL;
  }
  return rm;
}

static Str* map_table(RegionMapping *rm, const char *sequence_region,
                      Error *err)
{
  Str *result = NULL;
  int has_err = 0;
  error_check(err);
  assert(rm && sequence_region);
  lua_pushstring(rm->L, sequence_region);
  lua_gettable(rm->L, -2); /* get mapping[sequence_region] */
  /* make sure mapping[sequence_region] is defined */
  if (lua_isnil(rm->L, -1)) {
    error_set(err, "mapping[%s] is nil (defined in \"%s\")", sequence_region,
              str_get(rm->mapping_filename));
    has_err = -1;
  }
  /* make sure mapping[sequence_region] is a string */
  if (!has_err) {
    if (!(lua_isstring(rm->L, -1))) {
      error_set(err, "mapping[%s] is not a string (defined in \"%s\")",
                sequence_region, str_get(rm->mapping_filename));
      has_err = -1;
    }
  }
  if (!has_err)
    result = str_new_cstr(lua_tostring(rm->L, -1));
  lua_pop(rm->L, 1); /* pop result */
  return result;
}

static Str* map_function(RegionMapping *rm, const char *sequence_region,
                         Error *err)
{
  Str *result = NULL;
  int has_err = 0;
  error_check(err);
  assert(rm && sequence_region);
  lua_getglobal(rm->L, "mapping");
  lua_pushstring(rm->L, sequence_region);
  /* call function */
  if (lua_pcall(rm->L, 1, 1, 0)) {
    error_set(err, "running function 'mapping': %s", lua_tostring(rm->L, -1));
    has_err = -1;
  }
  /* make sure the result is a string */
  if (!has_err) {
    if (!lua_isstring(rm->L, -1)) {
      error_set(err, "function 'mapping' must return a string (defined in "
                "\"%s\")", str_get(rm->mapping_filename));
      has_err = -1;
    }
    if (!has_err)
      result = str_new_cstr(lua_tostring(rm->L, -1));
    lua_pop(rm->L, 1); /* pop result */
  }
  return result;
}

Str* regionmapping_map(RegionMapping *rm, const char *sequence_region,
                       Error *err)
{
  error_check(err);
  assert(rm && sequence_region);
  if (rm->is_table)
    return map_table(rm, sequence_region, err);
  else
    return map_function(rm, sequence_region, err);
}

void regionmapping_free(RegionMapping *rm)
{
  if (!rm) return;
  str_free(rm->mapping_filename);
  if (rm->L) lua_close(rm->L);
  free(rm);
}
