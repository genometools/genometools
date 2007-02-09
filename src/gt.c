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
  int rval;
  gtr = gtr_new();
  gtr_register_components(gtr);
  rval = gtr_parse(gtr, argc, argv);
  argc -= rval;
  argv += rval;
  rval = gtr_run(gtr, argc, argv);
  gtr_free(gtr);
  return rval;
}
