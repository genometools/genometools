/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <string.h>
#include "error.h"
#include "fileutils.h"
#include "gtdata.h"
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "splitter.h"
#include "xansi.h"

#define GTDATADIR "/gtdata"
#define UPDIR     "/.."

Str* gtdata_get_path(const char *prog)
{
  Str *path;
  assert(prog);
  path = file_find_in_path(prog);
  assert(path);
  str_append_cstr(path, GTDATADIR);
  if (file_exists(str_get(path)))
    return path;
  str_set_length(path, str_length(path) - strlen(GTDATADIR));
  str_append_cstr(path, UPDIR);
  str_append_cstr(path, GTDATADIR);
  if (file_exists(str_get(path)))
    return path;
  error("could not find gtdata/ directory");
  assert(0);
}

void gtdata_show_help(const char *progname, void *unused)
{
  Splitter *splitter;
  Str *doc_file;
  lua_State *L;
  char *prog;

  assert(progname);

  prog = xstrdup(progname); /* create modifiable copy for splitter */
  splitter = splitter_new();
  splitter_split(splitter, prog, strlen(prog), ' ');
  doc_file = gtdata_get_path(splitter_get_token(splitter, 0));
  assert(doc_file);
  str_append_cstr(doc_file, "/doc/");

  /* create Lua & push gtdata_doc_dir to Lua */
  L = luaL_newstate();
  luaL_openlibs(L);
  if (!L)
    error("out of memory (cannot create new lua state)");
  lua_pushstring(L, str_get(doc_file));
  lua_setglobal(L, "gtdata_doc_dir");

  /* finish creating doc_file */
  str_append_cstr(doc_file,
                  splitter_get_token(splitter, splitter_size(splitter) - 1));
  str_append_cstr(doc_file, ".lua");

  /* execute doc_file */
  if (luaL_loadfile(L, str_get(doc_file)) || lua_pcall(L, 0, 0, 0))
    error("cannot run doc file: %s", lua_tostring(L, -1));

  /* free */
  lua_close(L);
  str_free(doc_file);
  splitter_free(splitter);
  free(prog);
}
