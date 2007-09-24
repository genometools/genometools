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
#include "libgtmatch/esa-seqread.h"
#include "libgtmatch/esa-mmsearch-def.h"
#include "libgtmatch/format64.h"
#include "libgtmatch/verbose-def.h"

#include "libgtmatch/esa-maxpairs.pr"
#include "libgtmatch/test-maxpairs.pr"

typedef struct
{
  unsigned int userdefinedleastlength;
  unsigned long samples;
  bool scanfile;
  Str *indexname;
  StrArray *queryfiles;
} Maxpairsoptions;

static int simpleexactselfmatchoutput(/*@unused@*/ void *info,
                                      Seqpos len,
                                      Seqpos pos1,
                                      Seqpos pos2,
                                      /*@unused@*/ Env *env)
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

static int simpleexactquerymatchoutput(/*@unused@*/ void *info,
                                       unsigned long len,
                                       Seqpos dbstart,
                                       uint64_t unitnum,
                                       unsigned long querystart,
                                       /*@unused@*/ Env *env)
{
  printf("%lu " FormatSeqpos " " Formatuint64_t " %lu\n",
           len,
           PRINTSeqposcast(dbstart),
           PRINTuint64_tcast(unitnum),
           querystart);
  return 0;
}

static OPrval parse_options(Maxpairsoptions *maxpairsoptions,
                            int *parsed_args,
                            int argc,
                            const char **argv,
                            Env *env)
{
  OptionParser *op;
  Option *option, *queryoption, *scanoption, *sampleoption;
  OPrval oprval;

  env_error_check(env);
  op = option_parser_new("[options] -ii indexname",
                         "Perform Substring matches with or without query.",
                         env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  option = option_new_uint_min("l","Specify minimum length",
                               &maxpairsoptions->userdefinedleastlength,
                               (unsigned int) 20,
                               (unsigned int) 1,env);
  option_parser_add_option(op, option, env);

  sampleoption = option_new_ulong_min("samples","Specify number of samples",
                                 &maxpairsoptions->samples,
                                 (unsigned long) 0,
                                 (unsigned long) 1,
                                 env);
  option_parser_add_option(op, sampleoption, env);

  scanoption = option_new_bool("scan","scan index",
                               &maxpairsoptions->scanfile,
                               false,
                               env);
  option_parser_add_option(op, scanoption, env);

  option = option_new_string("ii",
                             "Specify input index",
                             maxpairsoptions->indexname, NULL, env);
  option_parser_add_option(op, option, env);
  option_is_mandatory(option);

  queryoption = option_new_filenamearray("q",
                             "Specify query files",
                             maxpairsoptions->queryfiles, env);
  option_parser_add_option(op, queryoption, env);

  oprval = option_parser_parse(op, parsed_args, argc, argv,
                               versionfunc, env);
  if (option_is_set(queryoption))
  {
    if (option_is_set(sampleoption))
    {
      env_error_set(env,"option -samples cannot be combined with option -q");
      return OPTIONPARSER_ERROR;
    }
    if (option_is_set(scanoption))
    {
      env_error_set(env,"option -scan cannot be combined with option -q");
      return OPTIONPARSER_ERROR;
    }
  }
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
  maxpairsoptions.queryfiles = strarray_new(env);
  oprval = parse_options(&maxpairsoptions,&parsed_args, argc, argv, env);
  if (oprval == OPTIONPARSER_OK)
  {
    Verboseinfo *verboseinfo = newverboseinfo(false,env);
    assert(parsed_args == argc);
    if (strarray_size(maxpairsoptions.queryfiles) == 0)
    {
      if (maxpairsoptions.samples == 0)
      {
        if (callenummaxpairs(maxpairsoptions.indexname,
                             maxpairsoptions.userdefinedleastlength,
                             maxpairsoptions.scanfile,
                             simpleexactselfmatchoutput,
                             NULL,
                             verboseinfo,
                             env) != 0)
        {
          haserr = true;
        }
      } else
      {
        if (testmaxpairs(maxpairsoptions.indexname,
                         maxpairsoptions.samples,
                         maxpairsoptions.userdefinedleastlength,
                         (Seqpos) (100 *
                                   maxpairsoptions.userdefinedleastlength),
                         verboseinfo,
                         env) != 0)
        {
          haserr = true;
        }
      }
    } else
    {
      if (callenumquerymatches(maxpairsoptions.indexname,
                               maxpairsoptions.queryfiles,
                               false,
                               maxpairsoptions.userdefinedleastlength,
                               simpleexactquerymatchoutput,
                               NULL,
                               verboseinfo,
                               env) != 0)
      {
        haserr = true;
      }
    }
    freeverboseinfo(&verboseinfo,env);
  }
  str_delete(maxpairsoptions.indexname,env);
  strarray_delete(maxpairsoptions.queryfiles,env);

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
