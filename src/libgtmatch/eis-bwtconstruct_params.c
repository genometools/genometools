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

#include "libgtcore/error.h"
#include "libgtcore/option.h"
#include "libgtmatch/eis-bwtseq.h"

#include "libgtmatch/eis-blockenc_params.h"
#include "libgtmatch/eis-bwtconstruct_params.h"

extern void
registerPackedIndexOptions(OptionParser *op, struct bwtOptions *paramOutput,
                           int defaultOptimizationFlags,
                           const Str *projectName)
{
  Option *option;

  registerBlockEncOptions(op, &paramOutput->final.seqParams.blockEnc);

  paramOutput->final.baseType = BWT_ON_BLOCK_ENC;

  paramOutput->final.featureToggles = BWTBaseFeatures;

  option = option_new_uint(
    "locfreq", "specify the locate frequency\n"
    "parameter i means that each i-th position of input string is stored\n"
    "0 => no locate information", &paramOutput->final.locateInterval, 16U);
  option_parser_add_option(op, option);

  option = option_new_bool(
    "locbitmap", "marked/unmarked positions for locate are stored as bitmaps\n"
    "this gives faster location of hits but increases the index by 1 bit per"
    " symbol", &paramOutput->useLocateBitmap, true);
  option_parser_add_option(op, option);
  paramOutput->useLocateBitmapOption = option;

  paramOutput->final.projectName = projectName;

  paramOutput->defaultOptimizationFlags = defaultOptimizationFlags;
}

static int
estimateBestLocateTypeFeature(const struct bwtOptions *paramOutput, Error *err)
{
  /* two cases: we store 1 bit per position or log(segmentlen) for
   * each marked position plus one to note the number of marked positions */
  unsigned segmentLen = estimateSegmentSize(&paramOutput->final.seqParams,
                                            paramOutput->final.baseType, err);
  if (segmentLen > (segmentLen + 1) * requiredUIntBits(segmentLen)
      / paramOutput->final.locateInterval)
    return BWTLocateCount;
  else
    return BWTLocateBitmap;
}

extern void
computePackedIndexDefaults(struct bwtOptions *paramOutput, int extraToggles,
                           Error *err)
{
  if (option_is_set(paramOutput->useLocateBitmapOption))
    paramOutput->final.featureToggles
      |= (paramOutput->useLocateBitmap?BWTLocateBitmap:BWTLocateCount);
  else
    paramOutput->final.featureToggles
      |= estimateBestLocateTypeFeature(paramOutput, err);
  paramOutput->final.featureToggles |= extraToggles;
  switch (paramOutput->final.baseType)
  {
  case BWT_ON_BLOCK_ENC:
    {
      int *EISFeatureSet = &paramOutput->final.seqParams.blockEnc.EISFeatureSet;
      *EISFeatureSet = EIS_FEATURE_NONE;
      if (paramOutput->defaultOptimizationFlags & BWTDEFOPT_LOW_RAM_OVERHEAD)
      {
        *EISFeatureSet &= ~EIS_FEATURE_REGION_SUMS;
      }
      if (paramOutput->defaultOptimizationFlags & BWTDEFOPT_FAST_RANK)
      {
        *EISFeatureSet |= EIS_FEATURE_REGION_SUMS;
      }
    }
    break;
  default:
    break;
  }
}
