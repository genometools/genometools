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
#include "core/assert_api.h"
#include "core/fileutils.h"
#include "core/gtdatapath.h"

#define GTDATADIR "/gtdata"
#define UPDIR     "/.."
static const char* GTDATA_DEFAULT_PATHS[]={ "/usr/share/genometools" GTDATADIR,
                                            NULL };

GtStr* gt_get_gtdata_path(const char *prog, GtError *err)
{
  GtStr *path;
  const char **defaultpath;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(prog);
  path = gt_str_new();
  had_err = gt_file_find_exec_in_path(path, prog, err);
  if (!had_err) {
    gt_assert(gt_str_length(path));
    gt_str_append_cstr(path, GTDATADIR);
    if (gt_file_exists(gt_str_get(path)))
      return path;
    gt_str_set_length(path, gt_str_length(path) - strlen(GTDATADIR));
    gt_str_append_cstr(path, UPDIR);
    gt_str_append_cstr(path, GTDATADIR);
    if (gt_file_exists(gt_str_get(path)))
      return path;
    for (defaultpath = GTDATA_DEFAULT_PATHS; *defaultpath; defaultpath++) {
      gt_str_reset(path);
      gt_str_append_cstr(path, *defaultpath);
      if (gt_file_exists(gt_str_get(path)))
        return path;
    }
    if (!gt_file_exists(gt_str_get(path))) {
      gt_error_set(err, "could not find gtdata/ directory");
      had_err = -1;
    }
  }
  if (had_err) {
    gt_str_delete(path);
    return NULL;
  }
  return path;
}
