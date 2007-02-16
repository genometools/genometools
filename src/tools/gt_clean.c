/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static void remove_pattern_in_current_dir(const char *pattern)
{
  char **files_to_remove;
  Str *path;
  glob_t g;

  path = str_new_cstr("./*");
  str_append_cstr(path, pattern);
  xglob(str_get(path), GLOB_NOCHECK, NULL, &g);

  /* remove found files */
  if (g.gl_pathc) {
    files_to_remove = g.gl_pathv;
    if (strcmp(*files_to_remove, str_get(path))) {
      while (*files_to_remove) {
        xunlink(*files_to_remove);
        files_to_remove++;
      }
    }
  }

  /* free */
  globfree(&g);
  str_free(path);
}

int gt_clean(int argc, char *argv[], Error *err)
{
  error_check(err);

  if (argc > 1) {
    fprintf(stderr, "Usage: %s\n", argv[0]);
    fprintf(stderr, "Remove all files in the current directory which are "
                    "automatically created by gt.\n");
    exit(EXIT_FAILURE); /* XXX */
  }

  /* remove GT_BIOSEQ_INDEX files */
  remove_pattern_in_current_dir(GT_BIOSEQ_INDEX);

  /* remove GT_BIOSEQ_RAW files */
  remove_pattern_in_current_dir(GT_BIOSEQ_RAW);

  return 0;
}
