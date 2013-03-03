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
#include "core/error.h"
#include "core/logger.h"
#include "core/minmax.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/versionfunc.h"
#include "match/eis-bwtseq.h"
#include "match/eis-bwtseq-param.h"
#include "core/encseq.h"
#include "match/enum-patt.h"
#include "match/esa-mmsearch.h"
#include "match/sarr-def.h"
#include "match/esa-map.h"
#include "tools/gt_packedindex_chk_search.h"
#include "match/sfx-apfxlen.h"

#define DEFAULT_PROGRESS_INTERVAL  100000UL

struct chkSearchOptions
{
  struct bwtOptions idx;
  long minPatLen, maxPatLen;
  unsigned long numOfSamples, progressInterval;
  int flags;
  bool verboseOutput;
};

static GtOPrval
parseChkBWTOptions(int *parsed_args, int argc, const char **argv,
                   struct chkSearchOptions *params, const GtStr *projectName,
                   GtError *err);

extern int
gt_packedindex_chk_search(int argc, const char *argv[], GtError *err)
{
  struct chkSearchOptions params;
  Suffixarray suffixarray;
  Enumpatterniterator *epi = NULL;
  bool saIsLoaded = false;
  BWTSeq *bwtSeq = NULL;
  GtStr *inputProject = NULL;
  int parsedArgs;
  bool had_err = false;
  BWTSeqExactMatchesIterator EMIter;
  bool EMIterInitialized = false;
  GtLogger *logger = NULL;
  inputProject = gt_str_new();

  do {
    gt_error_check(err);
    {
      bool exitNow = false;
      switch (parseChkBWTOptions(&parsedArgs, argc, argv, &params,
                                 inputProject, err))
      {
      case GT_OPTION_PARSER_OK:
        break;
      case GT_OPTION_PARSER_ERROR:
        had_err = true;
        exitNow = true;
        break;
      case GT_OPTION_PARSER_REQUESTS_EXIT:
        exitNow = true;
        break;
      }
      if (exitNow)
        break;
    }
    gt_str_set(inputProject, argv[parsedArgs]);

    logger = gt_logger_new(params.verboseOutput,
                           GT_LOGGER_DEFLT_PREFIX, stdout);

    bwtSeq = gt_availBWTSeq(&params.idx.final, logger, err);
    if ((had_err = bwtSeq == NULL))
      break;

    {
      enum verifyBWTSeqErrCode retval =
        gt_BWTSeqVerifyIntegrity(bwtSeq, gt_str_get(inputProject), params.flags,
                              params.progressInterval, stderr, logger, err);
      if ((had_err = (retval != VERIFY_BWTSEQ_NO_ERROR)))
      {
        fprintf(stderr, "index integrity check failed: %s\n",
                gt_error_get(err));
        gt_error_set(err, "aborted because of index integrity check fail");
        break;
      }
    }
    if (BWTSeqHasLocateInformation(bwtSeq))
    {
      if ((had_err = !gt_initEmptyEMIterator(&EMIter, bwtSeq)))
      {
        gt_error_set(err, "Cannot create matches iterator for sequence index.");
        break;
      }
      EMIterInitialized = true;
    }
    {
      unsigned long totalLen, dbstart;
      unsigned long trial, patternLen;

      if ((had_err =
           gt_mapsuffixarray(&suffixarray, SARR_SUFTAB | SARR_ESQTAB,
                             gt_str_get(inputProject), NULL, err) != 0))
      {
        gt_error_set(err, "Can't load suffix array project with"
                  " demand for encoded sequence and suffix table files\n");
        break;
      }
      totalLen = gt_encseq_total_length(suffixarray.encseq);
      saIsLoaded = true;
      if ((had_err = (params.minPatLen >= 0L && params.maxPatLen >= 0L
                      && params.minPatLen > params.maxPatLen)))
      {
        gt_error_set(err, "Invalid pattern lengths selected: min=%ld, max=%ld;"
                  " min <= max is required.", params.minPatLen,
                  params.maxPatLen);
        break;
      }
      if (params.minPatLen < 0 || params.maxPatLen < 0)
      {
        unsigned int numofchars
          = gt_alphabet_num_of_chars(
                               gt_encseq_alphabet(suffixarray.encseq));
        if (params.minPatLen < 0)
          params.minPatLen
            = gt_recommendedprefixlength(numofchars,
                                         totalLen,
                                         GT_RECOMMENDED_MULTIPLIER_DEFAULT,
                                         true);
        if (params.maxPatLen < 0)
          params.maxPatLen =
            MAX(params.minPatLen,
                125 * gt_recommendedprefixlength(numofchars,totalLen,
                                         GT_RECOMMENDED_MULTIPLIER_DEFAULT,
                                         true)/100);
        else
          params.maxPatLen = MAX(params.maxPatLen, params.minPatLen);
      }
      fprintf(stderr, "Using patterns of lengths %lu to %lu\n",
              params.minPatLen, params.maxPatLen);
      if ((had_err = totalLen + 1 != BWTSeqLength(bwtSeq)))
      {
        gt_error_set(err, "base suffix array and index have diferrent lengths!"
                          "%lu vs. %lu",  totalLen + 1,
                  BWTSeqLength(bwtSeq));
        break;
      }
      if ((had_err =
           (epi = gt_newenumpatterniterator(params.minPatLen, params.maxPatLen,
                                         suffixarray.encseq,
                                         err)) == NULL))
      {
        fputs("Creation of pattern iterator failed!\n", stderr);
        break;
      }
      for (trial = 0; !had_err && trial < params.numOfSamples; ++trial)
      {
        const GtUchar *pptr = gt_nextEnumpatterniterator(&patternLen, epi);
        GtMMsearchiterator *mmsi =
          gt_mmsearchiterator_new_complete_olain(suffixarray.encseq,
                                            suffixarray.suftab,
                                            0,  /* leftbound */
                                            totalLen, /* rightbound */
                                            0, /* offset */
                                            suffixarray.readmode,
                                            pptr,
                                            patternLen);
        if (BWTSeqHasLocateInformation(bwtSeq))
        {
          if ((had_err = !gt_reinitEMIterator(&EMIter, bwtSeq, pptr, patternLen,
                                           false)))
          {
            fputs("Internal error: failed to reinitialize pattern match"
                  " iterator", stderr);
            abort();
          }
          gt_assert(gt_EMINumMatchesTotal(&EMIter) ==
                    gt_BWTSeqMatchCount(bwtSeq, pptr, patternLen,
                                        false));
          gt_assert(gt_EMINumMatchesTotal(&EMIter)
                      == gt_mmsearchiterator_count(mmsi));
          while (gt_mmsearchiterator_next(&dbstart,mmsi))
          {
            unsigned long matchPos = 0;
            bool match = EMIGetNextMatch(&EMIter, &matchPos, bwtSeq);
            if ((had_err = !match))
            {
              gt_error_set(err,
                           "matches of packedindex expired before mmsearch!");
              break;
            }
            if ((had_err = matchPos != dbstart))
            {
              gt_error_set(err, "packedindex match doesn't equal mmsearch "
                           "match result!\n%lu vs. %lu\n",
                           matchPos, dbstart);
            }
          }
          if (!had_err)
          {
            unsigned long matchPos;
            bool trailingMatch = EMIGetNextMatch(&EMIter, &matchPos, bwtSeq);
            if ((had_err = trailingMatch))
            {
              gt_error_set(err, "matches of mmsearch expired before fmindex!");
              break;
            }
          }
        }
        else
        {
          unsigned long numFMIMatches = gt_BWTSeqMatchCount(bwtSeq, pptr,
                                                         patternLen,
                                                         false),
            numMMSearchMatches = gt_mmsearchiterator_count(mmsi);
          if ((had_err = numFMIMatches != numMMSearchMatches))
          {
            gt_error_set(err, "Number of matches not equal for suffix array ("
                              "%lu) and fmindex (%lu).\n",
                      numFMIMatches, numMMSearchMatches);
          }
        }
        gt_mmsearchiterator_delete(mmsi);
        mmsi = NULL;
        if (params.progressInterval && !((trial + 1) % params.progressInterval))
          putc('.', stderr);
      }
      if (params.progressInterval)
        putc('\n', stderr);
      fprintf(stderr, "Finished %lu of %lu matchings successfully.\n",
              trial, params.numOfSamples);
    }
  } while (0);
  if (EMIterInitialized) gt_destructEMIterator(&EMIter);
  if (saIsLoaded) gt_freesuffixarray(&suffixarray);
  gt_freeEnumpatterniterator(epi);
  if (bwtSeq) gt_deleteBWTSeq(bwtSeq);
  if (logger) gt_logger_delete(logger);
  if (inputProject) gt_str_delete(inputProject);
  return had_err?-1:0;
}

static GtOPrval
parseChkBWTOptions(int *parsed_args, int argc, const char **argv,
                   struct chkSearchOptions *params, const GtStr *projectName,
                   GtError *err)
{
  GtOptionParser *op;
  GtOPrval oprval;
  GtOption *option, *optionProgress;
  bool checkSuffixArrayValues, tryContextRetrieve, tryFullRegen;

  gt_error_check(err);
  op = gt_option_parser_new("indexname",
                         "Load (or build if necessary) BWT index for project"
                         " <indexname> and perform verification of search"
                         " results.");

  gt_registerPackedIndexOptions(op, &params->idx, BWTDEFOPT_MULTI_QUERY,
                             projectName);

  option = gt_option_new_long("minpatlen",
                           "minimum length of patterns searched for, -1 "
                           "implies automatic choice based on index "
                           "properties", &params->minPatLen, -1);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_long("maxpatlen",
                           "maximum length of patterns searched for, -1 "
                           "implies automatic choice based on index "
                           "properties", &params->maxPatLen, -1);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong("nsamples",
                            "number of sequences to search for",
                            &params->numOfSamples, 1000);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("chksfxarray",
                           "verify integrity of stored suffix array positions",
                           &checkSuffixArrayValues, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("full-lfmap",
                           "verify complete backwards regeneration of "
                           "original sequence", &tryFullRegen, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("chkcontext",
                           "verify integrity of regenerated sequence context",
                           &tryContextRetrieve, false);
  gt_option_parser_add_option(op, option);

  optionProgress = gt_option_new_ulong("ticks",
                                    "print dot after this many symbols "
                                    "tested okay", &params->progressInterval,
                                    DEFAULT_PROGRESS_INTERVAL);
  gt_option_parser_add_option(op, optionProgress);

  option = gt_option_new_bool("v",
                           "print verbose progress information",
                           &params->verboseOutput,
                           false);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_max_args(op, 1, 1);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);

  /* condense boolean options to flags field */
  params->flags = (checkSuffixArrayValues?VERIFY_BWTSEQ_SUFVAL:0)
    | (tryFullRegen?VERIFY_BWTSEQ_LFMAPWALK:0)
    | (tryContextRetrieve?VERIFY_BWTSEQ_CONTEXT:0);
  /* compute parameters currently not set from command-line or
   * determined indirectly */
  gt_computePackedIndexDefaults(&params->idx, BWTBaseFeatures);

  gt_option_parser_delete(op);

  return oprval;
}
