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

#include <string.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "core/assert_api.h"
#include "core/cstr_api.h"
#include "core/fileutils_api.h"
#include "core/basename_api.h"
#include "core/gtdatapath.h"
#include "core/ma.h"
#include "core/splitter.h"
#include "core/unused_api.h"
#include "extended/gtdatahelp.h"

int gt_gtdata_show_help(const char *progname, GT_UNUSED void *unused,
                        GtError *err)
{
  GtSplitter *splitter;
  GtStr *doc_file;
  lua_State *L = NULL;
  char *prog, *bn;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(progname);

  prog = gt_cstr_dup(progname); /* create modifiable copy for splitter */
  splitter = gt_splitter_new();
  gt_splitter_split(splitter, prog, strlen(prog), ' ');
  doc_file = gt_get_gtdata_path(gt_splitter_get_token(splitter, 0), err);
  if (!doc_file)
    had_err = -1;

  if (!had_err) {
    gt_str_append_cstr(doc_file, "/doc/");
    /* create Lua & push gtdata_doc_dir to Lua */
    L = luaL_newstate();
    if (!L) {
      gt_error_set(err, "out of memory (cannot create new Lua state)");
      had_err = -1;
    }
  }

  if (!had_err) {
    luaL_openlibs(L);
    lua_pushstring(L, gt_str_get(doc_file));
    lua_setglobal(L, "gtdata_doc_dir");
    /* finish creating doc_file */
    if (gt_splitter_size(splitter) == 1) {
      /* special case for `gt` */
      bn = gt_basename(progname);
      gt_str_append_cstr(doc_file, bn);
      gt_free(bn);
    }
    else {
      /* general case for the tools */
      gt_str_append_cstr(doc_file,
                      gt_splitter_get_token(splitter,
                                         gt_splitter_size(splitter) - 1));
    }
    gt_str_append_cstr(doc_file, ".lua");
    /* execute doc_file */
    if (luaL_loadfile(L, gt_str_get(doc_file)) || lua_pcall(L, 0, 0, 0)) {
      gt_error_set(err, "cannot run doc file: %s", lua_tostring(L, -1));
      had_err = -1;
    }
  }

  /* free */
  if (L) lua_close(L);
  gt_str_delete(doc_file);
  gt_splitter_delete(splitter);
  gt_free(prog);

  return had_err;
}
