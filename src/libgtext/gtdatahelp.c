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
#include <string.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "libgtcore/cstr.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/getbasename.h"
#include "libgtcore/gtdatapath.h"
#include "libgtcore/ma.h"
#include "libgtcore/splitter.h"
#include "libgtcore/unused.h"
#include "libgtext/gtdatahelp.h"

int gtdata_show_help(const char *progname, UNUSED void *unused, Error *err)
{
  Splitter *splitter;
  Str *doc_file;
  lua_State *L = NULL;
  char *prog, *bn;
  int had_err = 0;

  error_check(err);
  assert(progname);

  prog = cstr_dup(progname); /* create modifiable copy for splitter */
  splitter = splitter_new();
  splitter_split(splitter, prog, strlen(prog), ' ');
  doc_file = gtdata_get_path(splitter_get_token(splitter, 0), err);
  if (!doc_file)
    had_err = -1;

  if (!had_err) {
    str_append_cstr(doc_file, "/doc/");
    /* create Lua & push gtdata_doc_dir to Lua */
    L = luaL_newstate();
    if (!L) {
      error_set(err, "out of memory (cannot create new lua state)");
      had_err = -1;
    }
  }

  if (!had_err) {
    luaL_openlibs(L);
    lua_pushstring(L, str_get(doc_file));
    lua_setglobal(L, "gtdata_doc_dir");
    /* finish creating doc_file */
    if (splitter_size(splitter) == 1) {
      /* special case for `gt` */
      bn = getbasename(progname);
      str_append_cstr(doc_file, bn);
      ma_free(bn);
    }
    else {
      /* general case for the tools */
      str_append_cstr(doc_file,
                      splitter_get_token(splitter,
                                         splitter_size(splitter) - 1));
    }
    str_append_cstr(doc_file, ".lua");
    /* execute doc_file */
    if (luaL_loadfile(L, str_get(doc_file)) || lua_pcall(L, 0, 0, 0)) {
      error_set(err, "cannot run doc file: %s", lua_tostring(L, -1));
      had_err = -1;
    }
  }

  /* free */
  if (L) lua_close(L);
  str_delete(doc_file);
  splitter_delete(splitter);
  ma_free(prog);

  return had_err;
}
