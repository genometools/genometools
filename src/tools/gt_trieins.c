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

#include "libgtcore/versionfunc.h"
#include "libgtcore/option.h"
#include "libgtmatch/test-trieins.pr"

static OPrval parse_options(bool *onlyins,int *parsed_args,
                            int argc, const char **argv,Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;

  env_error_check(env);
  op = option_parser_new("[options] indexname",
                         "Perfrom trie insertions and check consistency.",
                         env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  option= option_new_bool("ins","perform only insertions",onlyins,false,env);
  option_parser_add_option(op, option, env);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, (unsigned int) 1,
                                            (unsigned int) 2, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_trieins(int argc, const char **argv, Env *env)
{
  Str *indexname;
  bool haserr = false;
  int parsed_args;
  bool onlyins = false;

  env_error_check(env);

  switch (parse_options(&onlyins,&parsed_args, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args == 1 || parsed_args == 2);

  indexname = str_new_cstr(argv[parsed_args],env);
  if (test_trieins(onlyins,indexname,env) != 0)
  {
    haserr = true;
  }
  str_delete(indexname);
  return haserr ? -1 : 0;
}
