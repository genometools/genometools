/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <inttypes.h>
#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/sarr-def.h"

#include "libgtmatch/sfx-map.pr"
#include "libgtmatch/esa-maxpairs.pr"

typedef struct
{
  unsigned int userdefinedleastlength;
  Str *indexname;
} Maxpairsoptions;

static int simpleexactselfmatchoutput(/*@unused@*/ void *info,
                                      Seqpos len,
                                      Seqpos pos1,
                                      Seqpos pos2)
{
  if (pos1 > pos2)
  {
    Seqpos tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  printf(FormatSeqpos " " FormatSeqpos " " FormatSeqpos "\n",
            PRINTSeqposcast(len),
            PRINTSeqposcast(pos1),
            PRINTSeqposcast(pos2));
  return 0;
}

static int callenummaxpairs(const Str *indexname,
                            uint32_t userdefinedleastlength,
                            Env *env)
{
  Suffixarray suffixarray;
  Seqpos totallength;
  bool haserr = false;

  if (streamsuffixarray(&suffixarray,
                        &totallength,
                        SARR_LCPTAB | SARR_SUFTAB | SARR_ESQTAB,
                        indexname,
                        false,
                        env) != 0)
  {
    haserr = true;
  }
  if (!haserr && enumeratemaxpairs(&suffixarray,
                                   userdefinedleastlength,
                                   simpleexactselfmatchoutput,
                                   NULL,
                                   env) != 0)
  {
    haserr = true;
  }
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}

static OPrval parse_options(Maxpairsoptions *maxpairsoptions,
                            int *parsed_args,
                            int argc, const char **argv,Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;

  env_error_check(env);
  op = option_parser_new("[-l minlength] -ii indexname",
                         "Enumerate maximal pairs of minimum length.",
                         env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  option = option_new_uint_min("l","Specify minimum length",
                               &maxpairsoptions->userdefinedleastlength,
                               (unsigned int) 20,
                               (unsigned int) 1,env);
  option_parser_add_option(op, option, env);

  option = option_new_string("ii",
                             "Specify input index",
                             maxpairsoptions->indexname, NULL, env);
  option_parser_add_option(op, option, env);
  option_is_mandatory(option);

  oprval = option_parser_parse(op, parsed_args, argc, argv,
                               versionfunc, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_maxpairs(int argc, const char **argv, Env *env)
{
  bool haserr = false;
  int parsed_args;
  Maxpairsoptions maxpairsoptions;
  OPrval oprval;

  env_error_check(env);

  maxpairsoptions.indexname = str_new(env);
  oprval = parse_options(&maxpairsoptions,&parsed_args, argc, argv, env);
  if (oprval == OPTIONPARSER_OK)
  {
    assert(parsed_args == argc);
    if (callenummaxpairs(maxpairsoptions.indexname,
                        (uint32_t) maxpairsoptions.userdefinedleastlength,
                        env) != 0)
    {
      haserr = true;
    }
  }
  str_delete(maxpairsoptions.indexname,env);
  if (oprval == OPTIONPARSER_REQUESTS_EXIT)
  {
    return 0;
  }
  if (oprval == OPTIONPARSER_ERROR)
  {
    return -1;
  }
  return haserr ? -1 : 0;
}
