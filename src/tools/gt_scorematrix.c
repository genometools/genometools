/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

int gt_scorematrix(int argc, char *argv[], Error *err)
{
  ScoreMatrix *sm;
  int has_err = 0;
  error_check(err);
  if (argc != 2) {
    fprintf(stderr, "Usage: %s scorematrix_filename\n", argv[0]);
    fprintf(stderr, "Parse the given protein score matrix and show it on "
                    "stdout.\n");
    exit(EXIT_FAILURE); /* XXX */
  }
  sm = scorematrix_read_protein(argv[1], err);
  if (!sm)
    has_err = -1;
  if (!has_err)
    scorematrix_show(sm, stdout);
  scorematrix_free(sm);
  return has_err;
}
