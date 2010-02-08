/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2009 Center for Bioinformatics, University of Hamburg

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

#include "gtr.h"

/* The GenomeTools (gt) genome analysis system */
int main(int argc, char *argv[])
{
  GtError *err;
  GtR *gtr;
  int rval;
  gt_lib_init();
  err = gt_error_new();
  gt_error_set_progname(err, argv[0]);
  if (!(gtr = gtr_new(err))) {
    fprintf(stderr, "%s: error: %s\n", gt_error_get_progname(err),
            gt_error_get(err));
    return EXIT_FAILURE;
  }
  gtr_register_components(gtr);
  switch (gtr_parse(gtr, &rval, argc, (const char**) argv, err)) {
    case GT_OPTION_PARSER_OK:
      argc -= rval;
      argv += rval;
      rval = gtr_run(gtr, argc, (const char**) argv, err);
      break;
    case GT_OPTION_PARSER_ERROR:
      rval = EXIT_FAILURE; /* user error */
      break;
    case GT_OPTION_PARSER_REQUESTS_EXIT:
      rval = EXIT_SUCCESS; /* everything went fine */
  }
  if (gt_error_is_set(err)) {
    fprintf(stderr, "%s: error: %s\n", gt_error_get_progname(err),
            gt_error_get(err));
    gt_assert(rval);
  }
  gtr_delete(gtr);
  gt_error_delete(err);
  if (gt_lib_clean())
    return GT_EXIT_PROGRAMMING_ERROR; /* programmer error */
  return rval;
}
