/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static OPrval parse_options(int *parsed_args, int argc, char **argv, Error *err)
{
  OptionParser *op;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("scorematrix_filename", "Parse the given protein "
                         "score matrix and show it on stdout.");
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, err);
  option_parser_free(op);
  return oprval;
}

int gt_scorematrix(int argc, char *argv[], Error *err)
{
  ScoreMatrix *sm;
  int parsed_args, has_err = 0;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args == 1);

  sm = scorematrix_read_protein(argv[1], err);
  if (!sm)
    has_err = -1;
  if (!has_err)
    scorematrix_show(sm, stdout);
  scorematrix_free(sm);

  return has_err;
}
