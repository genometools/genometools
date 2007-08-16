/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static OPrval parse_options(bool *usestream,bool *verbose,int *parsed_args, 
                            int argc, const char **argv,Env *env)
{
  OptionParser *op;
  Option *optionstream, *optionverbose;
  OPrval oprval;

  env_error_check(env);
  op = option_parser_new("[options] indexname", 
                         "Map or Stream <indexname> and check consistency.",
                         env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  optionstream = option_new_bool("stream","stream the index",
                                 usestream,false,env);
  option_parser_add_option(op, optionstream, env);
  optionverbose = option_new_bool("v","be verbose",verbose,false,env);
  option_parser_add_option(op, optionverbose, env);

  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 2, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_sfxmap(int argc, const char **argv, Env *env)
{
  Str *indexname;
  bool haserr = false;
  Suffixarray suffixarray;
  Seqpos totallength;
  int parsed_args;
  bool usestream = false,
       verbose = false;

  env_error_check(env);

  switch (parse_options(&usestream,&verbose,&parsed_args, argc, argv, env)) 
  {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args >= 1 && parsed_args <= 3);

  indexname = str_new_cstr(argv[parsed_args],env);
  if ((usestream ? streamsuffixarray : mapsuffixarray)(&suffixarray,
                                                       &totallength,
                                                       SARR_ALLTAB, 
                                                       indexname,
                                                       verbose,
                                                       env) != 0)
  {
    haserr = true;
  }
  str_delete(indexname,env);
  if (!haserr)
  {
    if(checkspecialranges(suffixarray.encseq,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    int readmode;
    
    for(readmode = 0; readmode < 4; readmode++)
    {
      if(isdnaalphabet(suffixarray.alpha) || 
         ((Readmode) readmode) == Forwardmode ||  
         ((Readmode) readmode) == Reversemode)
      /* if(((Readmode) readmode) == Forwardmode) */
      {
        if (testencodedsequence(suffixarray.filenametab,
                                suffixarray.encseq,
                                (Readmode) readmode,
                                getsymbolmapAlphabet(suffixarray.alpha),
                                env) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  if (!haserr)
  {
    if (verifymappedstr(&suffixarray,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && !usestream)
  {
    checkentiresuftab(suffixarray.encseq,
                      suffixarray.readmode,
                      getcharactersAlphabet(suffixarray.alpha),
                      suffixarray.suftab,
                      false, /* specialsareequal  */
                      true,  /* specialsareequalatdepth0 */
                      0,
                      env);
  }
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}
