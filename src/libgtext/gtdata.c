/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <string.h>
#include <libgtext/gtdata.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

#define GTDATADIR "/gtdata"
#define UPDIR     "/.."

Str* gtdata_get_path(const char *prog, Env *env)
{
  Str *path;
  int has_err = 0;
  env_error_check(env);
  assert(prog);
  path = str_new(env);
  has_err = file_find_in_path(path, prog, env);
  if (!has_err) {
    assert(str_length(path));
    str_append_cstr(path, GTDATADIR, env);
    if (file_exists(str_get(path)))
      return path;
    str_set_length(path, str_length(path) - strlen(GTDATADIR));
    str_append_cstr(path, UPDIR, env);
    str_append_cstr(path, GTDATADIR, env);
    if (!file_exists(str_get(path))) {
      env_error_set(env, "could not find gtdata/ directory");
      has_err = -1;
    }
  }
  if (has_err) {
    str_delete(path, env);
    return NULL;
  }
  return path;
}

int gtdata_show_help(const char *progname, void *unused, Env *env)
{
  Splitter *splitter;
  Str *doc_file;
  lua_State *L = NULL;
  char *prog;
  int has_err = 0;

  env_error_check(env);
  assert(progname);

  prog = cstr_dup(progname, env); /* create modifiable copy for splitter */
  splitter = splitter_new(env);
  splitter_split(splitter, prog, strlen(prog), ' ', env);
  doc_file = gtdata_get_path(splitter_get_token(splitter, 0), env);
  if (!doc_file)
    has_err = -1;

  if (!has_err) {
    str_append_cstr(doc_file, "/doc/", env);
    /* create Lua & push gtdata_doc_dir to Lua */
    L = luaL_newstate();
    if (!L) {
      env_error_set(env, "out of memory (cannot create new lua state)");
      has_err = -1;
    }
  }

  if (!has_err) {
    luaL_openlibs(L);
    lua_pushstring(L, str_get(doc_file));
    lua_setglobal(L, "gtdata_doc_dir");
    /* finish creating doc_file */
    str_append_cstr(doc_file,
                    splitter_get_token(splitter, splitter_size(splitter) - 1),
                    env);
    str_append_cstr(doc_file, ".lua", env);
    /* execute doc_file */
    if (luaL_loadfile(L, str_get(doc_file)) || lua_pcall(L, 0, 0, 0)) {
      env_error_set(env, "cannot run doc file: %s", lua_tostring(L, -1));
      has_err = -1;
    }
  }

  /* free */
  if (L) lua_close(L);
  str_delete(doc_file, env);
  splitter_delete(splitter, env);
  env_ma_free(prog, env);

  return has_err;
}
