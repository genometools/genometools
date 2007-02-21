/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <gt.h>
#include "gtr.h"

/* The GenomeTools (gt) genome analysis system */
int main(int argc, char *argv[])
{
  GTR *gtr;
  Error *err;
  int rval;
  gtr = gtr_new();
  err = error_new();
  gtr_register_components(gtr);
  switch (gtr_parse(gtr, &rval, argc, argv, err)) {
    case OPTIONPARSER_OK:
      argc -= rval;
      argv += rval;
      rval = gtr_run(gtr, argc, argv, err);
      break;
    case OPTIONPARSER_ERROR:
      rval = EXIT_FAILURE;
      break;
    case OPTIONPARSER_REQUESTS_EXIT:
      rval = EXIT_SUCCESS;
  }
  if (error_is_set(err)) {
    fprintf(stderr, "error: %s\n", error_get(err));
    assert(rval);
  }
  error_delete(err);
  gtr_free(gtr);
  return rval;
}
