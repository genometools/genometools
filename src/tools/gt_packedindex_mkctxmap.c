/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#include "libgtcore/error.h"
#include "libgtcore/option.h"
#include "libgtcore/str.h"
#include "libgtcore/versionfunc.h"

#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseq-construct.h"
#include "libgtmatch/eis-bwtseq-context-param.h"
#include "libgtmatch/eis-bwtseq-sass.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-map.pr"

#include "tools/gt_packedindex_mkctxmap.h"

struct mkCtxMapOptions
{
  int mapIntervalLog2;
  bool verboseOutput;
};

static OPrval
parseMkCtxMapOptions(int *parsed_args, int argc, const char **argv,
                     struct mkCtxMapOptions *params, Error *err);

extern int
gt_packedindex_mkctxmap(int argc, const char *argv[], Error *err)
{
  struct mkCtxMapOptions params;
  Str *projectName = NULL;
  Verboseinfo *verbosity = NULL;
  BWTSeq *bwtSeq = NULL;
  SASeqSrc *src;
  int parsedArgs;
  bool had_err = false;
  bool saInitialized = false, saiInitialized = false;
  Suffixarray sa;
  SuffixarrayFileInterface sai;
  projectName = str_new();

  do {
    error_check(err);
    {
      bool exitNow = false;
      switch (parseMkCtxMapOptions(&parsedArgs, argc, argv, &params, err))
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
    str_set(projectName, argv[parsedArgs]);
    verbosity = newverboseinfo(params.verboseOutput);
    /* try to find appropriate suffix source */
    {
      Seqpos len;
      if (streamsuffixarray(&sa, &len, SARR_SUFTAB, projectName, verbosity,
                            err))
      {
        error_unset(err);
        if (streamsuffixarray(&sa, &len, 0, projectName, verbosity, err))
        {
          had_err = true;
          break;
        }
        ++len;
        saInitialized = true;
        bwtSeq = loadBWTSeqForSA(projectName, BWT_ON_BLOCK_ENC,
                                 BWTDEFOPT_MULTI_QUERY, &sa, len, err);
        if (!(src = BWTSeqNewSASeqSrc(bwtSeq, NULL)))
        {
          error_set(err, "The project %s does not contain sufficient"
                    " information to regenerate the suffix array.",
                    str_get(projectName));
          had_err = true;
          break;
        }
      }
      else
      {
        ++len;
        saInitialized = true;
        initSuffixarrayFileInterface(&sai, len, &sa);
        src = SAI2SASS(&sai);
        saiInitialized = true;
      }
      {
        SeqDataReader readSfxIdx = SASSCreateReader(src, SFX_REQUEST_SUFTAB);
        BWTSeqContextRetriever *bwtSeqCR;
        BWTSeqContextRetrieverFactory *bwtSeqCRF
          = newBWTSeqContextRetrieverFactory(len, params.mapIntervalLog2);
        if (BWTSCRFReadAdvance(bwtSeqCRF, len, readSfxIdx, err)
            != len)
        {
          error_set(err, "Creation of context map unsuccessful: %s",
                    error_get(err));
          had_err = true;
          deleteBWTSeqContextRetrieverFactory(bwtSeqCRF);
          break;
        }
        bwtSeqCR = BWTSCRFGet(bwtSeqCRF, bwtSeq, projectName);
        deleteBWTSeqCR(bwtSeqCR);
        deleteBWTSeqContextRetrieverFactory(bwtSeqCRF);
      }
    }
  } while (0);
  if (bwtSeq)
  {
    SASSDelete(src);
    deleteBWTSeq(bwtSeq);
  }
  if (saiInitialized) destructSuffixarrayFileInterface(&sai);;
  if (saInitialized) freesuffixarray(&sa);
  if (verbosity) freeverboseinfo(&verbosity);
  if (projectName) str_delete(projectName);
  return had_err?-1:0;
}

static OPrval
parseMkCtxMapOptions(int *parsed_args, int argc, const char **argv,
                     struct mkCtxMapOptions *params, Error *err)
{
  OptionParser *op;
  OPrval oprval;
  Option *option;

  error_check(err);
  op = option_parser_new("indexname",
                         "Build BWT packedindex for project <indexname>.");
  registerCtxMapOptions(op, &params->mapIntervalLog2);

  option = option_new_bool("v",
                           "print verbose progress information",
                           &params->verboseOutput,
                           false);
  option_parser_add_option(op, option);

  option_parser_set_min_max_args(op, 1, 1);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);

  option_parser_delete(op);

  return oprval;
}
