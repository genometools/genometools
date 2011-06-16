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

#include "core/error.h"
#include "core/logger.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/versionfunc.h"
#include "match/eis-bwtseq.h"
#include "match/eis-bwtseq-construct.h"
#include "match/eis-bwtseq-context-param.h"
#include "match/eis-bwtseq-sass.h"

#include "match/sarr-def.h"
#include "match/esa-map.h"
#include "tools/gt_packedindex_mkctxmap.h"

struct mkCtxMapOptions
{
  int mapIntervalLog2;
  bool verboseOutput;
};

static GtOPrval
parseMkCtxMapOptions(int *parsed_args, int argc, const char **argv,
                     struct mkCtxMapOptions *params, GtError *err);

extern int
gt_packedindex_mkctxmap(int argc, const char *argv[], GtError *err)
{
  struct mkCtxMapOptions params;
  const char *projectName;
  GtLogger *logger = NULL;
  BWTSeq *bwtSeq = NULL;
  SASeqSrc *src;
  int parsedArgs;
  bool had_err = false;
  bool saInitialized = false, saiInitialized = false;
  Suffixarray sa;
  SuffixarrayFileInterface sai;

  do {
    gt_error_check(err);
    {
      bool exitNow = false;
      switch (parseMkCtxMapOptions(&parsedArgs, argc, argv, &params, err))
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
    projectName = argv[parsedArgs];
    logger = gt_logger_new(params.verboseOutput,
                           GT_LOGGER_DEFLT_PREFIX, stdout);
    /* try to find appropriate suffix source */
    {
      unsigned long len;
      if (streamsuffixarray(&sa, SARR_SUFTAB, projectName, logger, err))
      {
        gt_error_unset(err);
        if (streamsuffixarray(&sa, 0, projectName, logger, err))
        {
          had_err = true;
          break;
        }
        len = gt_encseq_total_length(sa.encseq) + 1;
        saInitialized = true;
        bwtSeq = gt_loadBWTSeqForSA(projectName, BWT_ON_BLOCK_ENC,
                                 BWTDEFOPT_MULTI_QUERY,
                                 gt_encseq_alphabet(sa.encseq), err);
        if (!(src = gt_BWTSeqNewSASeqSrc(bwtSeq, NULL)))
        {
          gt_error_set(err, "The project %s does not contain sufficient"
                       " information to regenerate the suffix array.",
                       projectName);
          had_err = true;
          break;
        }
      }
      else
      {
        len = gt_encseq_total_length(sa.encseq) + 1;
        saInitialized = true;
        gt_initSuffixarrayFileInterface(&sai, len, &sa);
        src = SAI2SASS(&sai);
        saiInitialized = true;
      }
      {
        SeqDataReader readSfxIdx = SASSCreateReader(src, SFX_REQUEST_SUFTAB);
        BWTSeqContextRetriever *bwtSeqCR;
        BWTSeqContextRetrieverFactory *bwtSeqCRF
          = gt_newBWTSeqContextRetrieverFactory(len, params.mapIntervalLog2);
        if (gt_BWTSCRFReadAdvance(bwtSeqCRF, len, readSfxIdx)
            != len)
        {
          gt_error_set(err, "Creation of context map unsuccessful: %s",
                    gt_error_get(err));
          had_err = true;
          gt_deleteBWTSeqContextRetrieverFactory(bwtSeqCRF);
          break;
        }
        bwtSeqCR = gt_BWTSCRFGet(bwtSeqCRF, bwtSeq, projectName);
        gt_deleteBWTSeqCR(bwtSeqCR);
        gt_deleteBWTSeqContextRetrieverFactory(bwtSeqCRF);
      }
    }
  } while (0);
  if (bwtSeq)
  {
    SASSDelete(src);
    gt_deleteBWTSeq(bwtSeq);
  }
  if (saiInitialized) gt_destructSuffixarrayFileInterface(&sai);;
  if (saInitialized) gt_freesuffixarray(&sa);
  if (logger) gt_logger_delete(logger);
  return had_err?-1:0;
}

static GtOPrval
parseMkCtxMapOptions(int *parsed_args, int argc, const char **argv,
                     struct mkCtxMapOptions *params, GtError *err)
{
  GtOptionParser *op;
  GtOPrval oprval;
  GtOption *option;

  gt_error_check(err);
  op = gt_option_parser_new("indexname",
                         "Build BWT packedindex for project <indexname>.");
  gt_registerCtxMapOptions(op, &params->mapIntervalLog2);

  option = gt_option_new_bool("v",
                           "print verbose progress information",
                           &params->verboseOutput,
                           false);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_max_args(op, 1, 1);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);

  gt_option_parser_delete(op);

  return oprval;
}
