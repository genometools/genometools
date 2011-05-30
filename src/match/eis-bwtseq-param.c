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

#include <limits.h>
#include "core/error.h"
#include "core/option_api.h"
#include "match/eis-bwtseq.h"
#include "match/eis-bwtseq-param.h"
#include "match/eis-bwtseq-context-param.h"
#include "match/eis-encidxseq-param.h"

void
gt_registerPackedIndexOptions(GtOptionParser *op,
                              struct bwtOptions *paramOutput,
                              int defaultOptimizationFlags,
                              const GtStr *projectName)
{
  GtOption *option;

  gt_registerEncIdxSeqOptions(op, &paramOutput->final.seqParams);

  paramOutput->final.featureToggles = BWTBaseFeatures;

  option = gt_option_new_uint(
    "locfreq", "specify the locate frequency\n"
    "parameter i means that each i-th position of input string is stored\n"
    "0 => no locate information", &paramOutput->final.locateInterval, 16U);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool(
    "locbitmap", "marked/unmarked positions for locate are stored as bitmaps\n"
    "this gives faster location of hits but increases the index by 1 bit per"
    " symbol", &paramOutput->useLocateBitmap, true);
  gt_option_parser_add_option(op, option);
  paramOutput->useLocateBitmapOption = option;

  option = gt_option_new_bool(
    "sprank", "build rank table for special symbols\n"
    "this produces an index which can be used to regenerate the "
    "original sequence but increases the memory used during index creation",
    &paramOutput->useSourceRank, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_int_min_max(
    "sprankilog", "specify the interval of rank sampling as log value\n"
    "parameter i means that each 2^i-th position of source is sampled for "
    "rank\nundefined => chooses default of log(log(sequence length))",
    &paramOutput->final.sourceRankInterval, -1, -1,
    sizeof (unsigned long) * CHAR_BIT - 1);
  gt_option_parser_add_option(op, option);

  gt_registerCtxMapOptions(op, &paramOutput->final.ctxMapILog);

  paramOutput->final.projectName = projectName;

  paramOutput->defaultOptimizationFlags = defaultOptimizationFlags;
}

static int
estimateBestLocateTypeFeature(const struct bwtOptions *paramOutput)
{
  if (!paramOutput->final.locateInterval)
    return BWTBaseFeatures;
  {
    /* two cases: we store 1 bit per position or log(segmentlen) for
     * each marked position plus one to note the number of marked positions */
    unsigned segmentLen = gt_estimateSegmentSize(&paramOutput->final.seqParams);
    if (segmentLen > (segmentLen + 1) * gt_requiredUIntBits(segmentLen)
                     / paramOutput->final.locateInterval)
      return BWTLocateCount;
    else
      return BWTLocateBitmap;
  }
}

void
gt_computePackedIndexDefaults(struct bwtOptions *paramOutput, int extraToggles)
{
  if (gt_option_is_set(paramOutput->useLocateBitmapOption))
    paramOutput->final.featureToggles
      |= (paramOutput->useLocateBitmap?BWTLocateBitmap:BWTLocateCount);
  else
    paramOutput->final.featureToggles
      |= estimateBestLocateTypeFeature(paramOutput);
  if (paramOutput->final.sourceRankInterval >= 0
      || paramOutput->useSourceRank)
    paramOutput->final.featureToggles |= BWTReversiblySorted;
  paramOutput->final.featureToggles |= extraToggles;
  paramOutput->final.seqParams.EISFeatureSet
    = gt_convertBWTOptFlags2EISFeatures(paramOutput->defaultOptimizationFlags);
}

int
gt_convertBWTOptFlags2EISFeatures(int BWTOptFlags)
{
  int EISFeatureSet = EIS_FEATURE_NONE;
  if (BWTOptFlags & BWTDEFOPT_LOW_RAM_OVERHEAD)
  {
    EISFeatureSet &= ~EIS_FEATURE_REGION_SUMS;
  }
  if (BWTOptFlags & BWTDEFOPT_FAST_RANK)
  {
    EISFeatureSet |= EIS_FEATURE_REGION_SUMS;
  }
  return EISFeatureSet;
}
