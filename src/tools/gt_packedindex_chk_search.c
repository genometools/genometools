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
#include "libgtmatch/eis-bwtconstruct_params.h"
#include "tools/gt_packedindex_chk_search.h"

#define DEFAULT_PROGRESS_INTERVAL  100000UL

struct chkSearchOptions
{
  struct bwtOptions idx;
  long minPatLen, maxPatLen;
  unsigned long numOfSamples, progressInterval;
  bool checkSuffixArrayValues;
};

static OPrval
parseChkBWTOptions(int *parsed_args, int argc, const char **argv,
                   struct chkSearchOptions *params, const Str *projectName,
                   Error *err);

extern int
gt_packedindex_chk_search(int argc, const char *argv[], Error *err)
{
  struct chkSearchOptions params;
  Suffixarray suffixarray;
  Enumpatterniterator *epi = NULL;
  bool saIsLoaded = false;
  struct BWTSeq *bwtSeq = NULL;
  Str *inputProject = NULL;
  int parsedArgs;
  bool had_err = false;

  inputProject = str_new();

  do {
    error_check(err);
    {
      bool exitNow = false;
      switch (parseChkBWTOptions(&parsedArgs, argc, argv, &params,
                                 inputProject, err))
      {
      case OPTIONPARSER_OK:
        break;
      case OPTIONPARSER_ERROR:
        had_err = true;
        exitNow = true;
        break;
      case OPTIONPARSER_REQUESTS_EXIT:
        exitNow = true;
        break;
      }
      if (exitNow)
        break;
    }
    str_set(inputProject, argv[parsedArgs]);
    {
      bwtSeq = availBWTSeq(&params.idx.final, err);
    }
    ensure(had_err, bwtSeq);
    if (had_err)
      break;
    if (params.checkSuffixArrayValues)
    {
      enum verifyBWTSeqErrCode retval =
        BWTSeqVerifyIntegrity(bwtSeq, inputProject, params.progressInterval,
                              stderr, err);
      if (retval != VERIFY_BWTSEQ_NO_ERROR)
      {
        fprintf(stderr, "index integrity check failed: %s\n",
                error_get(err));
        error_set(err, "aborted because of index integrity check fail");
        break;
      }
    }
    {
      Seqpos totalLen, dbstart;
      unsigned long trial, patternLen;
      ensure(had_err,
             mapsuffixarray(&suffixarray, &totalLen,
                            SARR_SUFTAB | SARR_ESQTAB,
                            inputProject, NULL, err) == 0);
      if (had_err)
      {
        error_set(err, "Can't load suffix array project with"
                       " demand for encoded sequence and suffix table files\n");
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
      ensure(had_err, totalLen + 1 == BWTSeqLength(bwtSeq));

      ensure(had_err,
             (epi = newenumpatterniterator(params.minPatLen, params.maxPatLen,
                                           suffixarray.encseq,
                                           getnumofcharsAlphabet(
                                              suffixarray.alpha),
                                           err)));
      if (had_err)
      {
        fputs("Creation of pattern iterator failed!\n", stderr);
        break;
      }
      for (trial = 0; !had_err && trial < params.numOfSamples; ++trial)
      {
        const Uchar *pptr = nextEnumpatterniterator(&patternLen, epi);
        MMsearchiterator *mmsi =
          newmmsearchiterator(suffixarray.encseq,
                              suffixarray.suftab,
                              0,  /* leftbound */
                              totalLen, /* rightbound */
                              0, /* offset */
                              suffixarray.readmode,
                              pptr,
                              patternLen);
        BWTSeqExactMatchesIterator *EMIter =
          newEMIterator(bwtSeq, pptr, patternLen, err);
        Seqpos numMatches = EMINumMatchesTotal(EMIter);
        ensure(had_err, EMIter);
        if (had_err)
          break;
        assert(numMatches == BWTSeqMatchCount(bwtSeq, pptr,
                                              patternLen, err));
        assert(EMINumMatchesTotal(EMIter) == countmmsearchiterator(mmsi));
        fprintf(stderr, "trial %lu, "FormatSeqpos" matches\n"
                "pattern: ", trial,
                numMatches);
        showsymbolstringgeneric(stderr, suffixarray.alpha, pptr, patternLen);
        fputs("\n", stderr);
        while (nextmmsearchiterator(&dbstart,mmsi))
        {
          struct MatchData *match =
            EMIGetNextMatch(EMIter, bwtSeq, err);
          ensure(had_err, match);
          if (had_err)
          {
            fputs("matches of fmindex expired before mmsearch!\n", stderr);
            break;
          }
          ensure(had_err, match->sfxArrayValue == dbstart);
          if (had_err)
          {
            fprintf(stderr, "fmindex match doesn't equal mmsearch match "
                    "result!\n"FormatSeqpos" vs. "FormatSeqpos"\n",
                    match->sfxArrayValue, dbstart);
            had_err = true;
          }
        }
        if (!had_err)
        {
          struct MatchData *trailingMatch =
            EMIGetNextMatch(EMIter, bwtSeq, err);
          ensure(had_err, !trailingMatch);
          if (had_err)
          {
            fputs("matches of mmsearch expired before fmindex!\n", stderr);
            break;
          }
        }
        deleteEMIterator(EMIter,err);
        freemmsearchiterator(&mmsi);
      }
      fprintf(stderr, "Finished %lu of %lu matchings successfully.\n",
              trial, params.numOfSamples);
    }
  } while (0);
  if (saIsLoaded) freesuffixarray(&suffixarray);
  if (epi) freeEnumpatterniterator(&epi);
  if (bwtSeq) deleteBWTSeq(bwtSeq, err);
  if (inputProject) str_delete(inputProject);
  return had_err?-1:0;
}

static OPrval
parseChkBWTOptions(int *parsed_args, int argc, const char **argv,
                   struct chkSearchOptions *params, const Str *projectName,
                   Error *err)
{
  OptionParser *op;
  OPrval oprval;
  Option *option;

  error_check(err);
  op = option_parser_new("indexname",
                         "Load (or build if necessary) BWT index for project"
                         " <indexname> and perform verification of search"
                         " results.");

  registerPackedIndexOptions(op, &params->idx, BWTDEFOPT_MULTI_QUERY,
                             projectName);

  option = option_new_long("minpatlen",
                           "minimum length of patterns searched for, -1 "
                           "implies automatic choice based on index "
                           "properties", &params->minPatLen, -1);
  option_parser_add_option(op, option);

  option = option_new_long("maxpatlen",
                           "maximum length of patterns searched for, -1 "
                           "implies automatic choice based on index "
                           "properties", &params->maxPatLen, -1);
  option_parser_add_option(op, option);

  option = option_new_ulong("nsamples",
                            "number of sequences to search for",
                            &params->numOfSamples, 1000);
  option_parser_add_option(op, option);

  option = option_new_bool("chksfxarray",
                           "verify integrity of stored suffix array positions",
                           &params->checkSuffixArrayValues, false);
  option_parser_add_option(op, option);

  option = option_new_ulong("ticks", "print dot after this many symbols"
                            " tested okay", &params->progressInterval,
                            DEFAULT_PROGRESS_INTERVAL);
  option_parser_add_option(op, option);

  oprval = option_parser_parse_min_max_args(op, parsed_args, argc,
                                            argv, versionfunc, 1, 1, err);
  /* compute parameters currently not set from command-line or
   * determined indirectly */
  computePackedIndexDefaults(&params->idx,
                             BWTBaseFeatures & ~BWTProperlySorted, err);

  option_parser_delete(op);

  return oprval;
}
