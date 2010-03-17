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

#include "core/log.h"
#include "match/sarr-def.h"
#include "match/esa-map.h"

#include "match/eis-bitpackseqpos.h"
#include "match/eis-bwtseq.h"
#include "match/eis-bwtseq-construct.h"
#include "match/eis-bwtseq-extinfo.h"
#include "match/eis-bwtseq-param.h"
#include "match/eis-bwtseq-priv.h"
#include "match/eis-encidxseq.h"
#include "match/eis-encidxseq-construct.h"

extern BWTSeq *
availBWTSeq(const struct bwtParam *params, GtLogger *verbosity,
            GtError *err)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  gt_assert(params && err);
  gt_error_check(err);
  if (streamsuffixarray(&suffixArray, SARR_SUFTAB | SARR_BWTTAB
                        | SARR_ESQTAB, params->projectName, verbosity, err))
  {
    gt_error_unset(err);
    if (streamsuffixarray(&suffixArray, SARR_SUFTAB | SARR_ESQTAB,
                          params->projectName, verbosity, err))
    {
      gt_error_unset(err);
      if (streamsuffixarray(&suffixArray, 0,
                            params->projectName, verbosity, err))
        return NULL;
    }
  }
  len = gt_encodedsequence_total_length(suffixArray.encseq) + 1;
  bwtSeq = availBWTSeqFromSA(params, &suffixArray, len, err);
  freesuffixarray(&suffixArray);
  return bwtSeq;
}

extern BWTSeq *
trSuftab2BWTSeq(const struct bwtParam *params, GtLogger *verbosity,
                GtError *err)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  gt_assert(params && err);
  gt_error_check(err);
  do
  {
    if (streamsuffixarray(&suffixArray,
                          SARR_SUFTAB | SARR_BWTTAB | SARR_ESQTAB,
                          params->projectName, verbosity, err))
    {
      gt_error_unset(err);
      if (streamsuffixarray(&suffixArray, SARR_SUFTAB | SARR_ESQTAB,
                            params->projectName, verbosity, err))
      {
        gt_error_set(err, "suffix array project %s does not hold required "
                     "suffix array (.suf) and encoded sequence (.esq) "
                     "information!", gt_str_get(params->projectName));
        break;
      }
    }
    len = gt_encodedsequence_total_length(suffixArray.encseq) + 1;
    bwtSeq = createBWTSeqFromSA(params, &suffixArray, len, err);
    freesuffixarray(&suffixArray);
  } while (0);
  return bwtSeq;
}

extern BWTSeq *
availBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                  Seqpos totalLen, GtError *err)
{
  BWTSeq *bwtSeq;
  gt_assert(sa && params && err);
  gt_error_check(err);
  /* try loading index */
  bwtSeq = loadBWTSeqForSA(params->projectName, params->seqParams.encType,
                           params->seqParams.EISFeatureSet,
                           sa, totalLen, err);
  /* if loading didn't work try on-demand creation */
  if (!bwtSeq)
  {
    gt_error_unset(err);
    bwtSeq = createBWTSeqFromSA(params, sa, totalLen, err);
  }
  else
  {
    fputs("Using pre-computed sequence index.\n", stderr);
  }
  return bwtSeq;
}

enum {
  GT_ALPHABETHANDLING_DEFAULT = 0,
  GT_ALPHABETHANDLING_W_RANK  = 1,
};

static const enum rangeSortMode GTAlphabetRangeSort[][2] =
{
  { SORTMODE_VALUE, SORTMODE_UNDEFINED },
  { SORTMODE_VALUE, SORTMODE_RANK }
};

extern BWTSeq *
loadBWTSeq(const GtStr *projectName, int BWTOptFlags, GtLogger *verbosity,
           GtError *err)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  gt_assert(projectName && err);
  gt_error_check(err);
  if (mapsuffixarray(&suffixArray, 0, projectName, verbosity, err))
    return NULL;
  len = gt_encodedsequence_total_length(suffixArray.encseq) + 1;
  bwtSeq = loadBWTSeqForSA(projectName, BWT_ON_BLOCK_ENC, BWTOptFlags,
                           &suffixArray, len, err);
  freesuffixarray(&suffixArray);
  return bwtSeq;
}

extern BWTSeq *
loadBWTSeqForSA(const GtStr *projectName, enum seqBaseEncoding encType,
                int BWTOptFlags, const Suffixarray *sa,
                Seqpos totalLen, GtError *err)
{
  struct BWTSeq *bwtSeq = NULL;
  EISeq *seqIdx = NULL;
  MRAEnc *alphabet = NULL;
  gt_assert(projectName && sa && err);
  alphabet = SANewMRAEnc(sa);
  if ((seqIdx = loadEncIdxSeqForSA(
         sa, totalLen, projectName, encType,
         convertBWTOptFlags2EISFeatures(BWTOptFlags), err)))
    bwtSeq = newBWTSeq(seqIdx, alphabet,
                       GTAlphabetRangeSort[GT_ALPHABETHANDLING_DEFAULT]);
  if (!bwtSeq)
  {
    MRAEncDelete(alphabet);
    if (seqIdx)
      deleteEncIdxSeq(seqIdx);
  }
  return bwtSeq;
}

extern BWTSeq *
createBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                   Seqpos totalLen, GtError *err)
{
  BWTSeq *bwtSeq = NULL;
  if (!sa->longest.defined)
  {
    gt_log_log("error: position of null-rotation/longest suffix not available"
            " for project %s\n", gt_str_get(params->projectName));
  }
  else
  {
    SuffixarrayFileInterface sai;
    initSuffixarrayFileInterface(&sai, totalLen, sa);
    bwtSeq = createBWTSeqFromSAI(params, &sai, err);
    destructSuffixarrayFileInterface(&sai);
  }
  return bwtSeq;
}

static inline void
buildSpRTable(const struct bwtParam *params,
              Seqpos totalLen,
              const GtEncodedsequence *encseq,
              GtReadmode readmode,
              SpecialsRankLookup **sprTable,
              const enum rangeSortMode **rangeSort)
{
  if (params->featureToggles & BWTReversiblySorted)
  {
    int sampleIntervalLog2 = params->sourceRankInterval;
    if (sampleIntervalLog2 == -1)
    {
      sampleIntervalLog2
        = gt_requiredUIntBits(requiredSeqposBits(totalLen));
    }
    *sprTable = newSpecialsRankLookup(encseq, readmode, sampleIntervalLog2);
  }
  *rangeSort = GTAlphabetRangeSort[sprTable?
                                   GT_ALPHABETHANDLING_W_RANK:
                                   GT_ALPHABETHANDLING_DEFAULT];
}

static BWTSeq *
createBWTSeqFromSASS(const struct bwtParam *params, SASeqSrc *src,
                     SpecialsRankLookup *sprTable,
                     const enum rangeSortMode *rangeSort,
                     GtError *err);

extern BWTSeq *
createBWTSeqFromSAI(const struct bwtParam *params,
                    SuffixarrayFileInterface *sai,
                    GtError *err)
{
  BWTSeq *bwtSeq;
  SpecialsRankLookup *sprTable = NULL;
  const enum rangeSortMode *rangeSort;
  gt_assert(sai && err && params);
  buildSpRTable(params, SAIGetLength(sai), SAIGetEncSeq(sai),
                SAIGetGtReadmode(sai), &sprTable, &rangeSort);
  bwtSeq = createBWTSeqFromSASS(params, SAI2SASS(sai), sprTable, rangeSort,
                                err);
  if (sprTable)
    deleteSpecialsRankLookup(sprTable);
  return bwtSeq;
}

extern BWTSeq *
createBWTSeqFromSfxI(const struct bwtParam *params, sfxInterface *sfxi,
                     GtError *err)
{
  BWTSeq *bwtSeq;
  SpecialsRankLookup *sprTable = NULL;
  const enum rangeSortMode *rangeSort;
  gt_assert(sfxi && params && err);
  buildSpRTable(params, SfxIGetLength(sfxi), SfxIGetEncSeq(sfxi),
                SfxIGetReadmode(sfxi), &sprTable, &rangeSort);
  bwtSeq = createBWTSeqFromSASS(params, SfxI2SASS(sfxi), sprTable, rangeSort,
                                err);
  if (sprTable)
    deleteSpecialsRankLookup(sprTable);
  return bwtSeq;
}

static BWTSeq *
createBWTSeqFromSASS(const struct bwtParam *params, SASeqSrc *src,
                     SpecialsRankLookup *sprTable,
                     const enum rangeSortMode *rangeSort,
                     GtError *err)
{
  EISeq *seqIdx = NULL;
  BWTSeq *bwtSeq = NULL;
  seqIdx = createBWTSeqGeneric(params, createEncIdxSeqGen, src,
                               rangeSort, sprTable, err);
  if (seqIdx)
  {
    MRAEnc *alphabet = SASSNewMRAEnc(src);
    bwtSeq = newBWTSeq(seqIdx, alphabet, rangeSort);
    if (!bwtSeq)
    {
      deleteEncIdxSeq(seqIdx);
      MRAEncDelete(alphabet);
    }
  }
  return bwtSeq;
}
