/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "gt.h"

static OPrval parse_options(Str *indexname,StrArray *indexnametab,
                            int *parsed_args, int argc,
                            const char **argv,Env *env)
{
  OptionParser *op;
  OPrval oprval;
  Option *option;

  env_error_check(env);
  op = option_parser_new("storeindex <mkvindex1> <mkvindex2> ...",
                         "Merge indexes into one index.",
                         env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  option = option_new_filenamearray("ii",
                                    "specify input index files (mandatory)",
                                    indexnametab,env);
  option_is_mandatory(option);
  option_parser_add_option(op, option, env);

  option = option_new_string("indexname",
                             "specify index to be created",
                             indexname, NULL, env);

  option_is_mandatory(option);
  option_parser_add_option(op, option, env);

  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_mergeesa(int argc, const char **argv, Env *env)
{
  Str *storeindex;
  StrArray *indexnametab;
  bool haserr = false;
  int parsed_args;

  env_error_check(env);

  storeindex = str_new(env);
  indexnametab = strarray_new(env);
  switch (parse_options(storeindex, indexnametab, &parsed_args,
                        argc, argv, env))
  {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
         haserr = true; break;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  if (!haserr)
  {
    unsigned long i;

    printf("# storeindex=%s\n",str_get(storeindex));
    for (i=0; i<strarray_size(indexnametab); i++)
    {
      printf("# input=%s\n",strarray_get(indexnametab,i));
    }
    if (performtheindexmerging(storeindex,
                              indexnametab,env) != 0)
    {
      haserr = true;
    }
  }
  str_delete(storeindex,env);
  strarray_delete(indexnametab,env);
  return haserr ? -1 : 0;
}
