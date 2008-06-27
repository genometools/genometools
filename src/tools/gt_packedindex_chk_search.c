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
#include "libgtcore/error.h"
#include "libgtcore/minmax.h"
#include "libgtcore/option.h"
#include "libgtcore/str.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseq-param.h"
#include "libgtmatch/encseq-def.h"
#include "libgtmatch/enum-patt-def.h"
#include "libgtmatch/esa-mmsearch-def.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-map.pr"
#include "libgtmatch/sfx-apfxlen.pr"
#include "libgtmatch/verbose-def.h"
#include "tools/gt_packedindex_chk_search.h"

#define DEFAULT_PROGRESS_INTERVAL  100000UL

struct chkSearchOptions
{
  struct bwtOptions idx;
  long minPatLen, maxPatLen;
  unsigned long numOfSamples, progressInterval;
  int flags;
  bool verboseOutput;
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
  BWTSeq *bwtSeq = NULL;
  Str *inputProject = NULL;
  int parsedArgs;
  bool had_err = false;
  BWTSeqExactMatchesIterator EMIter;
  bool EMIterInitialized = false;
  Verboseinfo *verbosity = NULL;
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

    verbosity = newverboseinfo(params.verboseOutput);

    bwtSeq = availBWTSeq(&params.idx.final, verbosity, err);
    if ((had_err = bwtSeq == NULL))
      break;

    {
      enum verifyBWTSeqErrCode retval =
        BWTSeqVerifyIntegrity(bwtSeq, inputProject, params.flags,
                              params.progressInterval, stderr, verbosity, err);
      if ((had_err = (retval != VERIFY_BWTSEQ_NO_ERROR)))
      {
        fprintf(stderr, "index integrity check failed: %s\n",
                error_get(err));
        error_set(err, "aborted because of index integrity check fail");
        break;
      }
    }
    if (BWTSeqHasLocateInformation(bwtSeq))
    {
      if ((had_err = !initEmptyEMIterator(&EMIter, bwtSeq)))
      {
        error_set(err, "Cannot create matches iterator for sequence index.");
        break;
      }
      EMIterInitialized = true;
    }
    {
      Seqpos totalLen, dbstart;
      unsigned long trial, patternLen;

      if ((had_err =
           mapsuffixarray(&suffixarray, &totalLen, SARR_SUFTAB | SARR_ESQTAB,
                          inputProject, NULL, err) != 0))
      {
        error_set(err, "Can't load suffix array project with"
                  " demand for encoded sequence and suffix table files\n");
        break;
      }
      saIsLoaded = true;
      if ((had_err = (params.minPatLen >= 0L && params.maxPatLen >= 0L
                      && params.minPatLen > params.maxPatLen)))
      {
        error_set(err, "Invalid pattern lengths selected: min=%ld, max=%ld;"
                  " min <= max is required.", params.minPatLen,
                  params.maxPatLen);
        break;
      }
      if (params.minPatLen < 0 || params.maxPatLen < 0)
      {
        unsigned int numofchars = getnumofcharsAlphabet(suffixarray.alpha);
        if (params.minPatLen < 0)
          params.minPatLen = recommendedprefixlength(numofchars, totalLen);
        if (params.maxPatLen < 0)
          params.maxPatLen =
            MAX(params.minPatLen,
                125 * recommendedprefixlength(numofchars, totalLen) / 100);
        else
          params.maxPatLen = MAX(params.maxPatLen, params.minPatLen);
      }
      fprintf(stderr, "Using patterns of lengths %lu to %lu\n",
              params.minPatLen, params.maxPatLen);
      if ((had_err = totalLen + 1 != BWTSeqLength(bwtSeq)))
      {
        error_set(err, "base suffix array and index have diferrent lengths!"
                  FormatSeqpos" vs. "FormatSeqpos,  totalLen + 1,
                  BWTSeqLength(bwtSeq));
        break;
      }
      if ((had_err =
           (epi = newenumpatterniterator(params.minPatLen, params.maxPatLen,
                                         suffixarray.encseq,
                                         getnumofcharsAlphabet(
                                           suffixarray.alpha),
                                         err)) == NULL))
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
        if (BWTSeqHasLocateInformation(bwtSeq))
        {
          Seqpos numMatches;
          if ((had_err = !reinitEMIterator(&EMIter, bwtSeq, pptr, patternLen,
                                           false)))
          {
            fputs("Internal error: failed to reinitialize pattern match"
                  " iterator", stderr);
            abort();
          }
          numMatches = EMINumMatchesTotal(&EMIter);
          assert(numMatches == BWTSeqMatchCount(bwtSeq, pptr, patternLen,
                                                false));
          assert(EMINumMatchesTotal(&EMIter) == countmmsearchiterator(mmsi));
/*        fprintf(stderr, "trial %lu, "FormatSeqpos" matches\n" */
/*                "pattern: ", trial, numMatches); */
/*        showsymbolstringgeneric(stderr, suffixarray.alpha, pptr, */
/*                                patternLen); */
/*        putc('\n', stderr); */
          while (nextmmsearchiterator(&dbstart,mmsi))
          {
            Seqpos matchPos = 0;
            bool match = EMIGetNextMatch(&EMIter, &matchPos, bwtSeq);
            if ((had_err = !match))
            {
              error_set(err, "matches of packedindex expired before mmsearch!");
              break;
            }
            if ((had_err = matchPos != dbstart))
            {
              error_set(err, "packedindex match doesn't equal mmsearch match "
                        "result!\n"FormatSeqpos" vs. "FormatSeqpos"\n",
                        matchPos, dbstart);
            }
          }
          if (!had_err)
          {
            Seqpos matchPos;
            bool trailingMatch = EMIGetNextMatch(&EMIter, &matchPos, bwtSeq);
            if ((had_err = trailingMatch))
            {
              error_set(err, "matches of mmsearch expired before fmindex!");
              break;
            }
          }
        }
        else
        {
          Seqpos numFMIMatches = BWTSeqMatchCount(bwtSeq, pptr, patternLen,
                                                  false),
            numMMSearchMatches = countmmsearchiterator(mmsi);
          if ((had_err = numFMIMatches != numMMSearchMatches))
          {
            error_set(err, "Number of matches not equal for suffix array ("
                      FormatSeqpos") and fmindex ("FormatSeqpos".\n",
                      numFMIMatches, numMMSearchMatches);
          }
        }
        freemmsearchiterator(&mmsi);
        if (params.progressInterval && !((trial + 1) % params.progressInterval))
          putc('.', stderr);
      }
      if (params.progressInterval)
        putc('\n', stderr);
      fprintf(stderr, "Finished %lu of %lu matchings successfully.\n",
              trial, params.numOfSamples);
    }
  } while (0);
  if (EMIterInitialized) destructEMIterator(&EMIter);
  if (saIsLoaded) freesuffixarray(&suffixarray);
  if (epi) freeEnumpatterniterator(&epi);
  if (bwtSeq) deleteBWTSeq(bwtSeq);
  if (verbosity) freeverboseinfo(&verbosity);
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
  Option *option, *optionProgress;
  bool checkSuffixArrayValues, tryContextRetrieve, tryFullRegen;

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
                           &checkSuffixArrayValues, false);
  option_parser_add_option(op, option);

  option = option_new_bool("full-lfmap",
                           "verify complete backwards regeneration of "
                           "original sequence", &tryFullRegen, false);
  option_parser_add_option(op, option);

  option = option_new_bool("chkcontext",
                           "verify integrity of regenerated sequence context",
                           &tryContextRetrieve, false);
  option_parser_add_option(op, option);

  optionProgress = option_new_ulong("ticks", "print dot after this many symbols"
                                    " tested okay", &params->progressInterval,
                                    DEFAULT_PROGRESS_INTERVAL);
  option_parser_add_option(op, optionProgress);

  option = option_new_bool("v",
                           "print verbose progress information",
                           &params->verboseOutput,
                           false);
  option_parser_add_option(op, option);

  option_parser_set_min_max_args(op, 1, 1);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);

  /* condense boolean options to flags field */
  params->flags = (checkSuffixArrayValues?VERIFY_BWTSEQ_SUFVAL:0)
    | (tryFullRegen?VERIFY_BWTSEQ_LFMAPWALK:0)
    | (tryContextRetrieve?VERIFY_BWTSEQ_CONTEXT:0);
  /* compute parameters currently not set from command-line or
   * determined indirectly */
  computePackedIndexDefaults(&params->idx, BWTBaseFeatures);

  option_parser_delete(op);

  return oprval;
}
