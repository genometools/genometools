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
#include "libgtcore/splitter.h"
#include "libgtcore/gtdatapath.h"

#define GTDATADIR "/gtdata"
#define UPDIR     "/.."

Str* gtdata_get_path(const char *prog, Error *err)
{
  Str *path;
  int had_err = 0;
  error_check(err);
  assert(prog);
  path = str_new();
  had_err = file_find_in_path(path, prog, err);
  if (!had_err) {
    assert(str_length(path));
    str_append_cstr(path, GTDATADIR);
    if (file_exists(str_get(path)))
      return path;
    str_set_length(path, str_length(path) - strlen(GTDATADIR));
    str_append_cstr(path, UPDIR);
    str_append_cstr(path, GTDATADIR);
    if (!file_exists(str_get(path))) {
      error_set(err, "could not find gtdata/ directory");
      had_err = -1;
    }
  }
  if (had_err) {
    str_delete(path);
    return NULL;
  }
  return path;
}
