/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg

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
  Error *err;
  GTR *gtr;
  int rval;
  allocators_init();
  err = error_new();
  error_set_progname(err, argv[0]);
  if (!(gtr = gtr_new(err))) {
    fprintf(stderr, "%s: error: %s\n", error_get_progname(err), error_get(err));
    return EXIT_FAILURE;
  }
  gtr_register_components(gtr);
  switch (gtr_parse(gtr, &rval, argc, (const char**) argv, err)) {
    case OPTIONPARSER_OK:
      argc -= rval;
      argv += rval;
      rval = gtr_run(gtr, argc, (const char**) argv, err);
      break;
    case OPTIONPARSER_ERROR:
      rval = 1; /* user error */
      break;
    case OPTIONPARSER_REQUESTS_EXIT:
      rval = 0; /* everything went fine */
  }
  if (error_is_set(err)) {
    fprintf(stderr, "%s: error: %s\n", error_get_progname(err), error_get(err));
    assert(rval);
  }
  gtr_delete(gtr);
  error_delete(err);
  if (allocators_clean())
    return 2; /* programmer error */
  return rval;
}
