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

#include <inttypes.h>
#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-seqread.h"

#include "libgtmatch/sfx-map.pr"
#include "libgtmatch/esa-maxpairs.pr"
#include "libgtmatch/test-maxpairs.pr"

typedef struct
{
  unsigned int userdefinedleastlength;
  unsigned long samples;
  bool mapped;
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
                            unsigned int userdefinedleastlength,
                            bool mapped,
                            Env *env)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  env_error_check(env);
  ssar = newSequentialsuffixarrayreaderfromfile(indexname,
                                                SARR_LCPTAB |
                                                SARR_SUFTAB |
                                                SARR_ESQTAB,
                                                mapped 
                                                  ? SEQ_mappedboth : SEQ_scan,
                                                env);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr && 
      enumeratemaxpairs(ssar,
                        getnumofcharsAlphabet(
                                 alphabetSequentialsuffixarrayreader(ssar)),
                        encseqSequentialsuffixarrayreader(ssar),
                        readmodeSequentialsuffixarrayreader(ssar),
                        userdefinedleastlength,
                        simpleexactselfmatchoutput,
                        NULL,
                        env) != 0)
  {
    haserr = true;
  }
  if(ssar != NULL)
  {
    freeSequentialsuffixarrayreader(&ssar,env);
  }
  return haserr ? -1 : 0;
}

static OPrval parse_options(Maxpairsoptions *maxpairsoptions,
                            int *parsed_args,
                            int argc,
                            const char **argv,
                            Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;

  env_error_check(env);
  op = option_parser_new("[options] -ii indexname",
                         "Enumerate maximal pairs of minimum length.",
                         env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  option = option_new_uint_min("l","Specify minimum length",
                               &maxpairsoptions->userdefinedleastlength,
                               (unsigned int) 20,
                               (unsigned int) 1,env);
  option_parser_add_option(op, option, env);

  option = option_new_ulong_min("samples","Specify number of samples",
                                 &maxpairsoptions->samples,
                                 (unsigned long) 0,
                                 (unsigned long) 1,
                                 env);
  option_parser_add_option(op, option, env);

  option = option_new_bool("map","map index",
                           &maxpairsoptions->mapped,
                           false,
                           env);
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
                         maxpairsoptions.userdefinedleastlength,
                         maxpairsoptions.mapped,
                         env) != 0)
    {
      haserr = true;
    }
    if (!haserr && maxpairsoptions.samples > 0)
    {
      if (testmaxpairs(maxpairsoptions.indexname,
                      maxpairsoptions.samples,
                      maxpairsoptions.userdefinedleastlength,
                      10 * maxpairsoptions.userdefinedleastlength,
                      env) != 0)
      {
        haserr = true;
      }
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
