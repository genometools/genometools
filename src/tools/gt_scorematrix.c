/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

int gt_scorematrix(int argc, char *argv[])
{
  ScoreMatrix *sm;
  Error *err = error_new();
  if (argc != 2) {
    fprintf(stderr, "Usage: %s scorematrix_filename\n", argv[0]);
    fprintf(stderr, "Parse the given protein score matrix and show it on "
                    "stdout.\n");
    return EXIT_FAILURE;
  }
  sm = scorematrix_read_protein(argv[1], err);
  if (sm)
    scorematrix_show(sm, stdout);
  scorematrix_free(sm);
  error_abort(err);
  error_free(err);
  return EXIT_SUCCESS;
}
