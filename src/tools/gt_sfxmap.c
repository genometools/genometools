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

#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/verbose-def.h"
#include "libgtmatch/esa-map.pr"
#include "libgtmatch/test-encseq.pr"
#include "libgtmatch/pos2seqnum.pr"
#include "libgtmatch/test-mappedstr.pr"
#include "libgtmatch/sfx-suftaborder.pr"
#include "libgtmatch/echoseq.pr"

typedef struct
{
  bool usestream, verbose;
  unsigned long trials;
} Sfxmapoptions;

static OPrval parse_options(Sfxmapoptions *sfxmapoptions,
                            int *parsed_args,
                            int argc, 
                            const char **argv,
                            Env *env)
{
  OptionParser *op;
  Option *optionstream, *optionverbose, *optiontrials;
  OPrval oprval;

  env_error_check(env);
  op = option_parser_new("[options] indexname",
                         "Map or Stream <indexname> and check consistency.",
                         env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  optionstream = option_new_bool("stream","stream the index",
                                 &sfxmapoptions->usestream,false,env);
  option_parser_add_option(op, optionstream, env);
  optionverbose = option_new_bool("v","be verbose",&sfxmapoptions->verbose,
                                  false,env);
  option_parser_add_option(op, optionverbose, env);

  optiontrials = option_new_ulong("trials","specify number of trials",
                                  &sfxmapoptions->trials,0,
                                  env);
  option_parser_add_option(op, optiontrials, env);

  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, (unsigned int) 1,
                                            (unsigned int) 2, env);
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
  Verboseinfo *verboseinfo;
  Sfxmapoptions sfxmapoptions;

  env_error_check(env);

  switch (parse_options(&sfxmapoptions,&parsed_args, argc, argv, env))
  {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args >= 1 && parsed_args <= 3);

  indexname = str_new_cstr(argv[parsed_args],env);
  verboseinfo = newverboseinfo(sfxmapoptions.verbose,env);
  if ((sfxmapoptions.usestream ? streamsuffixarray 
                               : mapsuffixarray)(&suffixarray,
                                                 &totallength,
                                                 SARR_ALLTAB,
                                                 indexname,
                                                 verboseinfo,
                                                 env) != 0)
  {
    haserr = true;
  }
  freeverboseinfo(&verboseinfo,env);
  str_delete(indexname,env);
  if (!haserr)
  {
    if (checkspecialrangesfast(suffixarray.encseq,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (checkmarkpos(suffixarray.encseq,suffixarray.numofdbsequences,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    int readmode;

    for (readmode = 0; readmode < 4; readmode++)
    {
      if (isdnaalphabet(suffixarray.alpha,env) ||
         ((Readmode) readmode) == Forwardmode ||
         ((Readmode) readmode) == Reversemode)
      {
        if (testencodedsequence(suffixarray.filenametab,
                                suffixarray.encseq,
                                (Readmode) readmode,
                                getsymbolmapAlphabet(suffixarray.alpha),
                                sfxmapoptions.trials,
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
  if (!haserr && !sfxmapoptions.usestream)
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
  if (!haserr)
  {
    checkalldescriptions(suffixarray.destab,suffixarray.destablength,
                         suffixarray.numofdbsequences,env);
  }
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}
