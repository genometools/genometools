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

#include <assert.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "libgtcore/bioseq.h"
#include "libgtext/regionmapping.h"

struct RegionMapping {
  Str *mapping_filename,
      *sequence_filename,
      *sequence_file, /* the (current) sequence file */
      *sequence_name; /* the (current) sequence name */
  lua_State *L;
  bool is_table;
  Bioseq *bioseq; /* the current bioseq */
};

RegionMapping* regionmapping_new_mapping(Str *mapping_filename, Env *env)
{
  RegionMapping *rm;
  int had_err = 0;
  env_error_check(env);
  assert(mapping_filename);
  /* alloc */
  rm = env_ma_calloc(env, 1, sizeof (RegionMapping));
  rm->mapping_filename = str_ref(mapping_filename);
  /* create new lua state (i.e., interpreter) */
  rm->L = luaL_newstate();
  if (!rm->L) {
    /* XXX: -> assert? */
    env_error_set(env, "out of memory (cannot create new lua state)");
    had_err = -1;
  }
  /* load the standard libs into the lua interpreter */
  if (!had_err)
    luaL_openlibs(rm->L);
  /* try to load & run mapping file */
  if (!had_err) {
    if (luaL_loadfile(rm->L, str_get(mapping_filename)) ||
        lua_pcall(rm->L, 0, 0, 0)) {
      env_error_set(env, "cannot run mapping file: %s",
                    lua_tostring(rm->L, -1));
      had_err = -1;
    }
  }
  /* make sure a global 'mapping' variable is defined */
  if (!had_err) {
    lua_getglobal(rm->L, "mapping");
    if (lua_isnil(rm->L, -1)) {
      env_error_set(env, "'mapping' is not defined in \"%s\"",
                str_get(mapping_filename));
      had_err = -1;
    }
  }
  /* make sure it is either a table or a function */
  if (!had_err) {
    if (!(lua_istable(rm->L, -1) || lua_isfunction(rm->L, -1))) {
      env_error_set(env, "'mapping' must be either a table or a function "
                    "(defined in \"%s\")", str_get(mapping_filename));
      had_err = -1;
    }
  }
  /* remember if it is a table or a function */
  if (!had_err) {
    if (lua_istable(rm->L, -1))
      rm->is_table = true;
    else {
      rm->is_table = false;
      lua_pop(rm->L, 1);
    }
  }
  /* return */
  if (had_err) {
    regionmapping_delete(rm, env);
    return NULL;
  }
  return rm;
}

RegionMapping* regionmapping_new_seqfile(Str *sequence_filename, Env *env)
{
  RegionMapping *rm;
  assert(sequence_filename && env);
  rm = env_ma_calloc(env, 1, sizeof (RegionMapping));
  rm->sequence_filename = str_ref(sequence_filename);
  return rm;
}

static Str* map_table(RegionMapping *rm, const char *sequence_region,
                      Env *env)
{
  Str *result = NULL;
  int had_err = 0;
  env_error_check(env);
  assert(rm && sequence_region);
  lua_pushstring(rm->L, sequence_region);
  lua_gettable(rm->L, -2); /* get mapping[sequence_region] */
  /* make sure mapping[sequence_region] is defined */
  if (lua_isnil(rm->L, -1)) {
    env_error_set(env, "mapping[%s] is nil (defined in \"%s\")",
                  sequence_region, str_get(rm->mapping_filename));
    had_err = -1;
  }
  /* make sure mapping[sequence_region] is a string */
  if (!had_err) {
    if (!(lua_isstring(rm->L, -1))) {
      env_error_set(env, "mapping[%s] is not a string (defined in \"%s\")",
                sequence_region, str_get(rm->mapping_filename));
      had_err = -1;
    }
  }
  if (!had_err)
    result = str_new_cstr(lua_tostring(rm->L, -1), env);
  lua_pop(rm->L, 1); /* pop result */
  return result;
}

static Str* map_function(RegionMapping *rm, const char *sequence_region,
                         Env *env)
{
  Str *result = NULL;
  int had_err = 0;
  env_error_check(env);
  assert(rm && sequence_region);
  lua_getglobal(rm->L, "mapping");
  lua_pushstring(rm->L, sequence_region);
  /* call function */
  if (lua_pcall(rm->L, 1, 1, 0)) {
    env_error_set(env, "running function 'mapping': %s",
                  lua_tostring(rm->L, -1));
    had_err = -1;
  }
  /* make sure the result is a string */
  if (!had_err) {
    if (!lua_isstring(rm->L, -1)) {
      env_error_set(env, "function 'mapping' must return a string (defined in "
                "\"%s\")", str_get(rm->mapping_filename));
      had_err = -1;
    }
    if (!had_err)
      result = str_new_cstr(lua_tostring(rm->L, -1), env);
    lua_pop(rm->L, 1); /* pop result */
  }
  return result;
}

static Str* regionmapping_map(RegionMapping *rm, const char *sequence_region,
                              Env *env)
{
  env_error_check(env);
  assert(rm && sequence_region);
  if (rm->sequence_filename)
    return str_ref(rm->sequence_filename);
  else if (rm->is_table)
    return map_table(rm, sequence_region, env);
  else
    return map_function(rm, sequence_region, env);
}

static int update_bioseq_if_necessary(RegionMapping *rm, Str *seqid, Env *env)
{
  int had_err = 0;
  env_error_check(env);
  assert(rm && seqid);
  if (!rm->sequence_file || str_cmp(rm->sequence_name, seqid)) {
    str_delete(rm->sequence_file, env);
    rm->sequence_file = regionmapping_map(rm, str_get(seqid), env);
    if (!rm->sequence_file)
      had_err = -1;
    else {
      if (!rm->sequence_name)
        rm->sequence_name = str_new(env);
      else
        str_reset(rm->sequence_name);
      str_append_str(rm->sequence_name, seqid, env);
      bioseq_delete(rm->bioseq, env);
      rm->bioseq = bioseq_new_str(rm->sequence_file, env);
      if (!rm->bioseq)
        had_err = -1;
    }
  }
  return had_err;
}

int regionmapping_get_raw_sequence(RegionMapping *rm, const char **raw,
                                   Str *seqid, Env *env)
{
  int had_err = 0;
  env_error_check(env);
  assert(rm && seqid);
  had_err = update_bioseq_if_necessary(rm, seqid, env);
  if (!had_err)
    *raw = bioseq_get_raw_sequence(rm->bioseq);
  return had_err;
}

int regionmapping_get_raw_sequence_length(RegionMapping *rm,
                                          unsigned long *length, Str *seqid,
                                          Env *env)
{
  int had_err = 0;
  env_error_check(env);
  assert(rm && seqid);
  had_err = update_bioseq_if_necessary(rm, seqid, env);
  if (!had_err)
    *length = bioseq_get_raw_sequence_length(rm->bioseq);
  return had_err;
}

void regionmapping_delete(RegionMapping *rm, Env *env)
{
  if (!rm) return;
  str_delete(rm->mapping_filename, env);
  str_delete(rm->sequence_filename, env);
  str_delete(rm->sequence_file, env);
  str_delete(rm->sequence_name, env);
  if (rm->L) lua_close(rm->L);
  bioseq_delete(rm->bioseq, env);
  env_ma_free(rm, env);
}
