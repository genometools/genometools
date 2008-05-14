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
#include "libgtcore/log.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-map.pr"

#include "libgtmatch/eis-bwtconstruct_params.h"
#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseqconstruct.h"
#include "libgtmatch/eis-bwtseqcreate.h"
#include "libgtmatch/eis-bwtseqpriv.h"
#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-encidxseqconstruct.h"

extern BWTSeq *
availBWTSeq(const struct bwtParam *params, Verboseinfo *verbosity, Error *err)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  assert(params && err);
  error_check(err);
  if (streamsuffixarray(&suffixArray, &len, SARR_SUFTAB | SARR_BWTTAB
                        | SARR_ESQTAB, params->projectName, verbosity, err))
  {
    error_unset(err);
    if (streamsuffixarray(&suffixArray, &len, SARR_SUFTAB | SARR_ESQTAB,
                          params->projectName, verbosity, err))
    {
      error_unset(err);
      if (streamsuffixarray(&suffixArray, &len, 0,
                            params->projectName, verbosity, err))
        return NULL;
    }
  }
  ++len;
  bwtSeq = availBWTSeqFromSA(params, &suffixArray, len, err);
  freesuffixarray(&suffixArray);
  return bwtSeq;
}

extern BWTSeq *
trSuftab2BWTSeq(const struct bwtParam *params, Verboseinfo *verbosity,
                Error *err)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  assert(params && err);
  error_check(err);
  do
  {
    if (streamsuffixarray(&suffixArray, &len,
                          SARR_SUFTAB | SARR_BWTTAB | SARR_ESQTAB,
                          params->projectName, verbosity, err))
    {
      error_unset(err);
      if (streamsuffixarray(&suffixArray, &len, SARR_SUFTAB | SARR_ESQTAB,
                            params->projectName, verbosity, err))
      {
        error_set(err, "suffix array project %s does not hold required suffix"
                  " array (.suf) and encoded sequence (.esq) information!",
                  str_get(params->projectName));
        break;
      }
    }
    ++len;
    bwtSeq = createBWTSeqFromSA(params, &suffixArray, len, err);
    freesuffixarray(&suffixArray);
  } while (0);
  return bwtSeq;
}

extern BWTSeq *
availBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                  Seqpos totalLen, Error *err)
{
  BWTSeq *bwtSeq;
  assert(sa && params && err);
  error_check(err);
  /* try loading index */
  bwtSeq = loadBWTSeqForSA(params->projectName, params->baseType,
                           params->seqParams.blockEnc.EISFeatureSet,
                           sa, totalLen, err);
  /* if loading didn't work try on-demand creation */
  if (!bwtSeq)
  {
    error_unset(err);
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
loadBWTSeq(const Str *projectName, int BWTOptFlags, Verboseinfo *verbosity,
           Error *err)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  assert(projectName && err);
  error_check(err);
  if (mapsuffixarray(&suffixArray, &len, 0, projectName, verbosity, err))
    return NULL;
  ++len;
  bwtSeq = loadBWTSeqForSA(projectName, BWT_ON_BLOCK_ENC, BWTOptFlags,
                           &suffixArray, len, err);
  freesuffixarray(&suffixArray);
  return bwtSeq;
}

extern BWTSeq *
loadBWTSeqForSA(const Str *projectName, enum seqBaseEncoding baseType,
                int BWTOptFlags, const Suffixarray *sa,
                Seqpos totalLen, Error *err)
{
  struct BWTSeq *bwtSeq = NULL;
  EISeq *seqIdx = NULL;
  MRAEnc *alphabet = NULL;
  assert(projectName && sa && err);
  alphabet = newMRAEncFromSA(sa);
  switch (baseType)
  {
  case BWT_ON_BLOCK_ENC:
    if ((seqIdx = loadBlockEncIdxSeqForSA(
           sa, totalLen, projectName,
           convertBWTOptFlags2EISFeatures(BWTOptFlags), err)))
      bwtSeq = newBWTSeq(seqIdx, alphabet,
                         GTAlphabetRangeSort[GT_ALPHABETHANDLING_DEFAULT]);
    break;
  default:
    error_set(err, "Illegal/unknown/unimplemented encoding requested!");
    break;
  }
  if (!bwtSeq && seqIdx)
    deleteEncIdxSeq(seqIdx);
  return bwtSeq;
}

extern BWTSeq *
createBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                   Seqpos totalLen, Error *err)
{
  BWTSeq *bwtSeq = NULL;
  if (!sa->longest.defined)
  {
    log_log("error: position of null-rotation/longest suffix not available"
            " for project %s\n", str_get(params->projectName));
  }
  else
  {
    SuffixarrayFileInterface sai;
    initSuffixarrayFileInterface(&sai, sa);
    bwtSeq = createBWTSeqFromSAI(params, &sai, totalLen, err);
    destructSuffixarrayFileInterface(&sai);
  }
  return bwtSeq;
}

extern BWTSeq *
createBWTSeqFromSAI(const struct bwtParam *params,
                    SuffixarrayFileInterface *sai,
                    Seqpos totalLen, Error *err)
{
  BWTSeq *bwtSeq = NULL;
  MRAEnc *alphabet = NULL;
  SeqDataReader readSfxIdx = { NULL, NULL };
  assert(sai && err && params);
  alphabet = newMRAEncFromSAI(sai);
  if (params->locateInterval
      && !SDRIsValid(readSfxIdx
                               = SAIMakeReader(sai, SFX_REQUEST_SUFTAB)))
  {
    error_set(err, "error: locate sampling requested but not available"
              " for project %s\n", str_get(params->projectName));
  }
  else
  {
    EISeq *seqIdx = NULL;
    switch (params->baseType)
    {
      RandomSeqAccessor origSeqRead;
    case BWT_ON_BLOCK_ENC:
      origSeqRead.accessFunc = SAIGetOrigSeqSym;
      origSeqRead.state = sai;
      seqIdx =
        createBWTSeqGeneric(
          params, (indexCreateFunc)newBlockEncIdxSeqFromSAI, sai, totalLen,
          alphabet, NULL, GTAlphabetRangeSort[GT_ALPHABETHANDLING_DEFAULT],
          origSeqRead, readSfxIdx, NULL, reportSAILongest, sai, err);
      break;
    default:
      error_set(err, "Illegal/unknown/unimplemented encoding requested!");
      break;
    }
    if (seqIdx)
      bwtSeq = newBWTSeq(seqIdx, alphabet,
                         GTAlphabetRangeSort[GT_ALPHABETHANDLING_DEFAULT]);
    if (!bwtSeq)
    {
      if (seqIdx)
        deleteEncIdxSeq(seqIdx);
      if (alphabet)
        MRAEncDelete(alphabet);
    }
  }
  return bwtSeq;
}

extern BWTSeq *
createBWTSeqFromSfxI(const struct bwtParam *params, sfxInterface *sfxi,
                     Seqpos totalLen, Error *err)
{
  EISeq *seqIdx = NULL;
  BWTSeq *bwtSeq = NULL;
  MRAEnc *alphabet = NULL;
  SeqDataReader readSfxIdx = { NULL, NULL };
  SpecialsRankTable *sprTable = NULL;
  const enum rangeSortMode *rangeSort;
  assert(sfxi && params && err);
  if (params->locateInterval)
  {
    if (!SDRIsValid(readSfxIdx
                              = SfxIRegisterReader(sfxi,
                                                   SFX_REQUEST_SUFTAB)))
      return NULL;
  }
  if (params->featureToggles & BWTReversiblySorted)
  {
    int sampleIntervalLog2 = params->sourceRankInterval;
    if (sampleIntervalLog2 == -1)
    {
      sampleIntervalLog2
        = requiredUIntBits(requiredSeqposBits(totalLen));
    }
    sprTable = newSpecialsRankTable(SfxIGetEncSeq(sfxi),
                                    SfxIGetReadMode(sfxi),
                                    sampleIntervalLog2);
  }
  rangeSort = GTAlphabetRangeSort[sprTable?
                                  GT_ALPHABETHANDLING_W_RANK:
                                  GT_ALPHABETHANDLING_DEFAULT];
  alphabet = newMRAEncFromSfxI(sfxi);
  {
    RandomSeqAccessor origSeqAccess = { SfxIGetOrigSeq, sfxi };
    seqIdx= createBWTSeqGeneric(
      params, (indexCreateFunc)newBlockEncIdxSeqFromSfxI, sfxi, totalLen,
      alphabet, getSfxISeqStats(sfxi), rangeSort,
      origSeqAccess, readSfxIdx, sprTable, (reportLongest)getSfxILongestPos,
      sfxi, err);
  }
  if (seqIdx)
  {
    bwtSeq = newBWTSeq(seqIdx, alphabet, rangeSort);
  }
  if (!bwtSeq && seqIdx)
    deleteEncIdxSeq(seqIdx);
  if (sprTable)
    deleteSpecialsRankTable(sprTable);
  return bwtSeq;
}
