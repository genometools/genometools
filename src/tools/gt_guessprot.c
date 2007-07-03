/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  /* XXX: add one liner to describe this tool */
  op = option_parser_new("filenames", 
                         "guess if sequence in filenames is protein or DNA.",
                          env);
  oprval = option_parser_parse_min_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_guessprot(int argc, const char **argv, Env *env)
{
  int i, parsed_args, retval;
  StrArray *filenametab;

  env_error_check(env);

  switch (parse_options(&parsed_args, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  filenametab = strarray_new(env);
  for(i=parsed_args; i < argc; i++)
  {
    strarray_add_cstr(filenametab,argv[i],env);
  }
  retval = guessifproteinsequencestream(filenametab,env);
  if(retval < 0)
  {
    return -1;
  }
  if(retval == 1)
  {
    exit(1); /* XXX */
  } else
  {
    return 0;
  }
}
