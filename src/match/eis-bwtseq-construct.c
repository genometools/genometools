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

#include "core/encseq_metadata.h"
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

BWTSeq *
gt_availBWTSeq(const struct bwtParam *params, GtLogger *verbosity,
            GtError *err)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  unsigned long len;
  gt_assert(params && err);
  gt_error_check(err);
  if (streamsuffixarray(&suffixArray, SARR_SUFTAB | SARR_BWTTAB
                        | SARR_ESQTAB, gt_str_get(params->projectName),
                        verbosity, err))
  {
    gt_error_unset(err);
    if (streamsuffixarray(&suffixArray, SARR_SUFTAB | SARR_ESQTAB,
                          gt_str_get(params->projectName), verbosity, err))
    {
      GtEncseqMetadata *emd = NULL;
      gt_error_unset(err);
      emd = gt_encseq_metadata_new(gt_str_get(params->projectName), err);
      if (emd == NULL)
        return NULL;
      len = gt_encseq_metadata_total_length(emd);
      gt_encseq_metadata_delete(emd);
    } else {
      len = gt_encseq_total_length(suffixArray.encseq) + 1;
    }
  } else {
    len = gt_encseq_total_length(suffixArray.encseq) + 1;
  }
  bwtSeq = gt_availBWTSeqFromSA(params, &suffixArray, len, err);
  gt_freesuffixarray(&suffixArray);
  return bwtSeq;
}

BWTSeq *
gt_trSuftab2BWTSeq(const struct bwtParam *params, GtLogger *verbosity,
                GtError *err)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  unsigned long len;
  gt_assert(params && err);
  gt_error_check(err);
  do
  {
    if (streamsuffixarray(&suffixArray,
                          SARR_SUFTAB | SARR_BWTTAB | SARR_ESQTAB,
                          gt_str_get(params->projectName), verbosity, err))
    {
      gt_error_unset(err);
      if (streamsuffixarray(&suffixArray, SARR_SUFTAB | SARR_ESQTAB,
                            gt_str_get(params->projectName), verbosity, err))
      {
        gt_error_set(err, "suffix array project %s does not hold required "
                     "suffix array (.suf) and encoded sequence (.esq) "
                     "information!", gt_str_get(params->projectName));
        break;
      }
    }
    len = gt_encseq_total_length(suffixArray.encseq) + 1;
    bwtSeq = gt_createBWTSeqFromSA(params, &suffixArray, len, err);
    gt_freesuffixarray(&suffixArray);
  } while (0);
  return bwtSeq;
}

BWTSeq *
gt_availBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                  unsigned long totalLen, GtError *err)
{
  BWTSeq *bwtSeq;
  gt_assert(sa && params && err);
  gt_error_check(err);
  /* try loading index */
  bwtSeq = gt_loadBWTSeqForSA(gt_str_get(params->projectName),
                              params->seqParams.encType,
                              params->seqParams.EISFeatureSet,
                              gt_encseq_alphabet(sa->encseq), err);
  /* if loading didn't work try on-demand creation */
  if (!bwtSeq)
  {
    gt_error_unset(err);
    bwtSeq = gt_createBWTSeqFromSA(params, sa, totalLen, err);
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

BWTSeq *
gt_loadBWTSeq(const char *projectName, int BWTOptFlags,
              GT_UNUSED GtLogger *verbosity, GtError *err)
{
  struct BWTSeq *bwtSeq = NULL;
  GtEncseq *encseq;
  GtEncseqLoader *el;
  gt_assert(projectName && err);
  gt_error_check(err);

  el = gt_encseq_loader_new();
  gt_encseq_loader_do_not_require_sds_tab(el);
  gt_encseq_loader_do_not_require_des_tab(el);
  gt_encseq_loader_do_not_require_ssp_tab(el);
  encseq = gt_encseq_loader_load(el, projectName, err);
  gt_encseq_loader_delete(el);

  if (encseq == NULL)
    return NULL;
  bwtSeq = gt_loadBWTSeqForSA(projectName, BWT_ON_BLOCK_ENC, BWTOptFlags,
                           gt_encseq_alphabet(encseq), err);
  gt_encseq_delete(encseq);
  return bwtSeq;
}

BWTSeq *
gt_loadBWTSeqForSA(const char *projectName, enum seqBaseEncoding encType,
                   int BWTOptFlags, const GtAlphabet *gtalphabet,
                   GtError *err)
{
  struct BWTSeq *bwtSeq = NULL;
  EISeq *seqIdx = NULL;
  MRAEnc *alphabet = NULL;
  gt_assert(projectName && gtalphabet && err);
  alphabet = gt_SANewMRAEnc(gtalphabet);
  seqIdx = gt_loadEncIdxSeqForSA(gtalphabet, projectName, encType,
                                 gt_convertBWTOptFlags2EISFeatures(BWTOptFlags),
                                 err);
  if (seqIdx != NULL)
  {
    bwtSeq = gt_newBWTSeq(seqIdx, alphabet,
                          GTAlphabetRangeSort[GT_ALPHABETHANDLING_DEFAULT]);
  }
  if (!bwtSeq)
  {
    gt_MRAEncDelete(alphabet);
    if (seqIdx)
    {
      gt_deleteEncIdxSeq(seqIdx);
    }
  }
  return bwtSeq;
}

BWTSeq *
gt_createBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                   unsigned long totalLen, GtError *err)
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
    gt_initSuffixarrayFileInterface(&sai, totalLen, sa);
    bwtSeq = gt_createBWTSeqFromSAI(params, &sai, err);
    gt_destructSuffixarrayFileInterface(&sai);
  }
  return bwtSeq;
}

static inline void
buildSpRTable(const struct bwtParam *params,
              unsigned long totalLen,
              const GtEncseq *encseq,
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
        = gt_requiredUIntBits(requiredUlongBits(totalLen));
    }
    *sprTable = gt_newSpecialsRankLookup(encseq, readmode, sampleIntervalLog2);
  }
  *rangeSort = GTAlphabetRangeSort[sprTable?
                                   GT_ALPHABETHANDLING_W_RANK:
                                   GT_ALPHABETHANDLING_DEFAULT];
}

static BWTSeq *
gt_createBWTSeqFromSASS(const struct bwtParam *params, SASeqSrc *src,
                     SpecialsRankLookup *sprTable,
                     const enum rangeSortMode *rangeSort,
                     GtError *err);

BWTSeq *
gt_createBWTSeqFromSAI(const struct bwtParam *params,
                    SuffixarrayFileInterface *sai,
                    GtError *err)
{
  BWTSeq *bwtSeq;
  SpecialsRankLookup *sprTable = NULL;
  const enum rangeSortMode *rangeSort;
  gt_assert(sai && err && params);
  buildSpRTable(params, SAIGetLength(sai), SAIGetEncSeq(sai),
                SAIGetGtReadmode(sai), &sprTable, &rangeSort);
  bwtSeq = gt_createBWTSeqFromSASS(params, SAI2SASS(sai), sprTable, rangeSort,
                                err);
  if (sprTable)
    gt_deleteSpecialsRankLookup(sprTable);
  return bwtSeq;
}

BWTSeq *
gt_createBWTSeqFromSfxI(const struct bwtParam *params, sfxInterface *sfxi,
                     GtError *err)
{
  BWTSeq *bwtSeq;
  SpecialsRankLookup *sprTable = NULL;
  const enum rangeSortMode *rangeSort;
  gt_assert(sfxi && params && err);
  buildSpRTable(params, gt_SfxIGetLength(sfxi), gt_SfxIGetEncSeq(sfxi),
                gt_SfxIGetReadmode(sfxi), &sprTable, &rangeSort);
  bwtSeq = gt_createBWTSeqFromSASS(params, gt_SfxI2SASS(sfxi), sprTable,
                                   rangeSort,
                                   err);
  if (sprTable)
    gt_deleteSpecialsRankLookup(sprTable);
  return bwtSeq;
}

static BWTSeq *
gt_createBWTSeqFromSASS(const struct bwtParam *params, SASeqSrc *src,
                     SpecialsRankLookup *sprTable,
                     const enum rangeSortMode *rangeSort,
                     GtError *err)
{
  EISeq *seqIdx = NULL;
  BWTSeq *bwtSeq = NULL;
  seqIdx = gt_createBWTSeqGeneric(params, gt_createEncIdxSeqGen, src,
                               rangeSort, sprTable, err);
  if (seqIdx)
  {
    MRAEnc *alphabet = SASSNewMRAEnc(src);
    bwtSeq = gt_newBWTSeq(seqIdx, alphabet, rangeSort);
    if (!bwtSeq)
    {
      gt_deleteEncIdxSeq(seqIdx);
      gt_MRAEncDelete(alphabet);
    }
  }
  return bwtSeq;
}
