/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#include <stdio.h>
#include <string.h>

#include "gt_packedindex_chk_search.h"
#include "libgtcore/ensure.h"
#include "libgtcore/env.h"
#include "libgtcore/minmax.h"
#include "libgtcore/option.h"
#include "libgtcore/str.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/encseq-def.h"
#include "libgtmatch/enum-patt-def.h"
#include "libgtmatch/esa-mmsearch-def.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-map.pr"
#include "libgtmatch/sfx-apfxlen.pr"
#include "gt_packedindex_bwtconstruct_params.h"

/***************************************************************************
 * routines for matching tests
 ***************************************************************************/

struct chkSearchOptions
{
  struct bwtOptions idx;
  long minPatLen, maxPatLen;
  unsigned locateInterval;
  unsigned long numOfSamples;
};

static OPrval
parseChkBWTOptions(int *parsed_args, int argc, const char **argv,
                   struct chkSearchOptions *param, const Str *projectName,
                   Env *env);

extern int
gt_packedindex_chk_search(int argc, const char *argv[], Env *env)
{
  struct chkSearchOptions params;
  Suffixarray suffixarray;
  bool saIsLoaded = false;
  struct BWTSeq *bwtSeq = NULL;
  Str *inputProject = NULL;
  Enumpatterniterator *epi = NULL;
  MMsearchiterator *mmsi = NULL;
  struct BWTSeqExactMatchesIterator *EMIter = NULL;
  int parsedArgs;
  bool had_err = false;

  inputProject = str_new(env);

  do {
    env_error_check(env);
    switch (parseChkBWTOptions(&parsedArgs, argc, argv, &params,
                               inputProject, env))
    {
    case OPTIONPARSER_OK:
      break;
    case OPTIONPARSER_ERROR:
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      return 0;
    }
    str_set(inputProject, argv[parsedArgs], env);
    {
      struct bwtParam bwtParams;
      memcpy(&bwtParams.seqParams, &params.idx,
             sizeof (union bwtSeqParam));
      bwtParams.baseType = BWT_ON_BLOCK_ENC;
      bwtParams.bwtFeatureToggles = BWTLocateBitmap;
      bwtParams.locateInterval = params.locateInterval;
      bwtParams.projectName = inputProject;
      bwtSeq = availBWTSeq(&bwtParams, env);
    }
    ensure(had_err, bwtSeq);
    if (had_err)
      break;
    {
      const Uchar *pptr;
      Seqpos totalLen, dbstart;
      unsigned long trial, patternLen;
      ensure(had_err,
             mapsuffixarray(&suffixarray, &totalLen,
                            SARR_SUFTAB | SARR_ESQTAB,
                            inputProject, NULL, env));
      if (had_err)
      {
        env_error_set(env, "Can't load suffix array project with"
                      " demand for BWT and suffix table files\n");
        break;
      }
      saIsLoaded = true;
      ensure(had_err, params.minPatLen < 0L || params.maxPatLen < 0L
             || params.minPatLen < params.maxPatLen);
      if (had_err)
        break;
      if (params.minPatLen < 0 || params.maxPatLen < 0)
      {
        unsigned int numofchars = getnumofcharsAlphabet(suffixarray.alpha);
        if (params.minPatLen < 0)
          params.minPatLen = recommendedprefixlength(numofchars, totalLen);
        if (params.maxPatLen < 0)
          params.maxPatLen =
            MAX(params.minPatLen,
                125 * recommendedprefixlength(numofchars, totalLen) / 100);
      }
      fprintf(stderr, "Using patterns of lengths %lu to %lu\n",
              params.minPatLen, params.maxPatLen);
      ensure(had_err, totalLen = BWTSeqLength(bwtSeq));

      ensure(had_err,
             (epi = newenumpatterniterator(params.minPatLen, params.maxPatLen,
                                           suffixarray.encseq, env)));
      if (had_err)
      {
        fputs("Creation of pattern iterator failed!\n", stderr);
        break;
      }
      for (trial = 0; !had_err && trial < params.numOfSamples; ++trial)
      {
        pptr = nextEnumpatterniterator(&patternLen, epi);
        mmsi = newmmsearchiterator(suffixarray.encseq,
                                   suffixarray.suftab,
                                   0,  /* leftbound */
                                   totalLen, /* rightbound */
                                   0, /* offset */
                                   suffixarray.readmode,
                                   pptr,
                                   patternLen,
                                   env);

        EMIter = newEMIterator(bwtSeq, pptr, patternLen, env);
        ensure(had_err, EMIter);
        if (had_err)
          break;
        assert(EMINumMatchesTotal(EMIter)
               == BWTSeqMatchCount(bwtSeq, pptr, patternLen, env));
        assert(EMINumMatchesTotal(EMIter) == countmmsearchiterator(mmsi));
        while (nextmmsearchiterator(&dbstart,mmsi))
        {
          struct MatchData *match =
            EMIGetNextMatch(EMIter, bwtSeq, env);
          ensure(had_err, match);
          if (had_err)
          {
            fputs("matches of fmindex expired before mmsearch!\n", stderr);
            break;
          }
          ensure(had_err, match->sfxArrayValue == dbstart);
          if (had_err)
          {
            fputs("fmindex match doesn't equal mmsearch match result!\n",
                  stderr);
            break;
          }
        }
        {
          struct MatchData *trailingMatch =
            EMIGetNextMatch(EMIter, bwtSeq, env);
          ensure(had_err, !trailingMatch);
          if (had_err)
          {
            fputs("matches of mmsearch expired before fmindex!\n", stderr);
            break;
          }
        }
      }
      fprintf(stderr, "Finished %lu matchings successfully.\n",
              params.numOfSamples);
    }
  } while (0);
  if (epi) freeEnumpatterniterator(&epi,env);
  if (mmsi) freemmsearchiterator(&mmsi,env);
  if (EMIter) deleteEMIterator(EMIter,env);
  if (bwtSeq) deleteBWTSeq(bwtSeq, env);
  if (inputProject) str_delete(inputProject, env);
  return had_err?-1:0;
}

static OPrval
parseChkBWTOptions(int *parsed_args, int argc, const char **argv,
                   struct chkSearchOptions *param, const Str *projectName,
                   Env *env)
{
  OptionParser *op;
  OPrval oprval;
  Option *option;

  env_error_check(env);
  op = option_parser_new("indexname",
                         "Load (or build if necessary) BWT index for project"
                         " <indexname> and perform verification of search"
                         " results.", env);

  registerPackedIndexOptions(op, &param->idx, BWTDEFOPT_CONSTRUCTION,
                             projectName, env);

  option = option_new_long("minpatlen",
                           "minimum length of patterns searched for, -1 "
                           "implies automatic choice based on index "
                           "properties", &param->minPatLen, -1, env);
  option_parser_add_option(op, option, env);

  option = option_new_long("maxpatlen",
                           "maximum length of patterns searched for, -1 "
                           "implies automatic choice based on index "
                           "properties", &param->maxPatLen, -1, env);
  option_parser_add_option(op, option, env);

  option = option_new_ulong("nsamples",
                            "number of sequences to search for",
                            &param->numOfSamples, 1000, env);
  option_parser_add_option(op, option, env);

  oprval = option_parser_parse_min_max_args(op, parsed_args, argc,
                                            argv, versionfunc, 1, 1, env);
  option_parser_delete(op, env);

  /* compute parameters currently not set from command-line or
   * determined indirectly */
  computePackedIndexDefaults(&param->idx, env);

  return oprval;
}
