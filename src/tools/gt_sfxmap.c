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
  op = option_parser_new("indexname", "ADD ONE LINER.", env);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_sfxmap(int argc, const char **argv, Env *env)
{
  const char *indexname;
  bool haserr = false;
  Suffixarray suffixarray;
  int parsed_args;

  env_error_check(env);

  switch (parse_options(&parsed_args, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args == 1);

  indexname = argv[parsed_args];

  if (mapsuffixarray(&suffixarray,true,true,indexname,env) != 0)
  {
    haserr = true;
  }

  if (!haserr)
  {
    if (testencodedsequence((const char **) suffixarray.filenametab,
                           suffixarray.numoffiles,
                           suffixarray.encseq,
                           getsymbolmapAlphabet(suffixarray.alpha),
                           env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (verifymappedstr(&suffixarray,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    checkentiresuftab(suffixarray.encseq,
                      getcharactersAlphabet(suffixarray.alpha),
                      suffixarray.suftab,
                      false, /* specialsareequal  */
                      true,  /* specialsareequalatdepth0 */
                      getencseqtotallength(suffixarray.encseq),
                      0,
                      env);
  }
  if (!haserr)
  {
    freesuffixarray(&suffixarray,env);
  }
  return haserr ? -1 : 0;
}
