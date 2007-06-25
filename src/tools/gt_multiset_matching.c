/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] multiset_string text",
                         "Match multiset defined by multiset_string against "
                         "text.",env);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 2, 2, env);
  option_parser_delete(op, env);
  return oprval;
}

static void show_match(unsigned long pos, void *data)
{
  printf("%lu\n", pos + 1);
}

int gt_multiset_matching(int argc, const char **argv, Env *env)
{
  int parsed_args;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* matching */
  assert(parsed_args + 1 < argc);
  multiset_matching((unsigned char*) argv[parsed_args],
                    strlen(argv[parsed_args]),
                    (unsigned char*) argv[parsed_args+1],
                    strlen(argv[parsed_args+1]), NULL, show_match);

  return 0;
}
