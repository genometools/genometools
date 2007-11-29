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
#include "libgtmatch/stamp.h"
#include "libgtmatch/enum-patt-def.h"
#include "libgtmatch/esa-mmsearch-def.h"

#include "libgtmatch/esa-map.pr"

typedef struct
{
  unsigned long minpatternlen, maxpatternlen, numofsamples;
  bool showpatt;
  Str *indexname;
} Pmatchoptions;

static int callpatternmatcher(const Pmatchoptions *pmopt,Env *env)
{
  Suffixarray suffixarray;
  Seqpos totallength;
  bool haserr = false;
  const Uchar *pptr;
  unsigned long patternlen;
  MMsearchiterator *mmsi;

  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     SARR_SUFTAB | SARR_ESQTAB,
                     pmopt->indexname,
                     NULL,
                     env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    unsigned long trial;
    Seqpos dbstart;
    Enumpatterniterator *epi;

    epi = newenumpatterniterator(pmopt->minpatternlen,
                                 pmopt->maxpatternlen,
                                 suffixarray.encseq,
                                 getnumofcharsAlphabet(suffixarray.alpha),
                                 env);
    for (trial = 0; trial < pmopt->numofsamples; trial++)
    {
      pptr = nextEnumpatterniterator(&patternlen,epi);
      if (pmopt->showpatt)
      {
        showsymbolstring(suffixarray.alpha,pptr,patternlen);
        printf("\n");
      }
      mmsi = newmmsearchiterator(suffixarray.encseq,
                                 suffixarray.suftab,
                                 0,  /* leftbound */
                                 totallength, /* rightbound */
                                 0, /* offset */
                                 suffixarray.readmode,
                                 pptr,
                                 patternlen,
                                 env);
      while (nextmmsearchiterator(&dbstart,mmsi))
      {
        /* Nothing */;
      }
      freemmsearchiterator(&mmsi,env);
    }
    if (pmopt->showpatt)
    {
      showPatterndistribution(epi);
    }
    freeEnumpatterniterator(&epi,env);
  }
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}

static OPrval parse_options(Pmatchoptions *pmopt,
                            int *parsed_args,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;

  error_check(err);
  op = option_parser_new("[options] -ii indexname",
                         "Perform pattern matches.");
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  option = option_new_ulong("minpl","Specify minimum length of pattern",
                           &pmopt->minpatternlen,
                           (unsigned long) 20);
  option_parser_add_option(op, option);
  option = option_new_ulong("maxpl","Specify maximum length of pattern",
                            &pmopt->maxpatternlen,
                            (unsigned long) 30);
  option_parser_add_option(op, option);

  option = option_new_ulong("samples","Specify number of samples",
                            &pmopt->numofsamples,
                           (unsigned long) 100000);
  option_parser_add_option(op, option);

  option = option_new_bool("s","Show generated pattern",
                            &pmopt->showpatt,
                            false);
  option_parser_add_option(op, option);

  option = option_new_string("ii",
                             "Specify input index",
                             pmopt->indexname, NULL);
  option_parser_add_option(op, option);
  option_is_mandatory(option);

  oprval = option_parser_parse(op, parsed_args, argc, argv,
                               versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

int gt_patternmatch(int argc, const char **argv, Env *env)
{
  bool haserr = false;
  int parsed_args;
  Pmatchoptions pmopt;
  OPrval oprval;

  env_error_check(env);

  pmopt.indexname = str_new();
  oprval = parse_options(&pmopt,&parsed_args, argc, argv, env_error(env));
  if (oprval == OPTIONPARSER_OK)
  {
    assert(parsed_args == argc);
    if (callpatternmatcher(&pmopt,env) != 0)
    {
      haserr = true;
    }
  }
  str_delete(pmopt.indexname);
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
