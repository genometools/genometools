/*
  Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#include <stdlib.h>

#include "libgtcore/chardef.h"
#include "libgtmatch/eis-bitpackseqpos.h"
#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseqcreate.h"
#include "libgtmatch/eis-bwtseqpriv.h"
#include "libgtmatch/eis-mrangealphabet.h"

/**
 * \file eis-bwtseqcreate.c generic methods for bwt index creation
 */

enum {
  LOCATE_INFO_IN_INDEX_HEADERID = 1111,
};

enum {
  LOCATE_HEADER_SIZE = sizeof (struct locateHeader),
};

struct locateHeaderWriteInfo
{
  unsigned locateInterval;
  reportLongest getLongest;
  void *getLongestState;
  int featureToggles;
};

static int
writeLocateInfoHeader(FILE *fp, void *cbData)
{
  struct locateHeader headerData;
  DefinedSeqpos longest;
  const struct locateHeaderWriteInfo *headerSrc = cbData;
  assert(cbData);
  headerData.locateInterval = headerSrc->locateInterval;
  longest = headerSrc->getLongest(headerSrc->getLongestState);
  if (!longest.defined)
  {
    fputs("Invalid index construction: position of suffix 0 unknown!\n",
          stderr);
    abort();
  }
  headerData.longest = longest.valueseqpos;
  headerData.featureToggles = headerSrc->featureToggles;
  return fwrite(&headerData, sizeof (headerData), 1, fp);
}

/**
 * @return 0 on error
 */
extern int
readLocateInfoHeader(EISeq *seqIdx, struct locateHeader *headerData)
{
  FILE *fp;
  assert(seqIdx && headerData);
  if (!(fp = EISSeekToHeader(seqIdx, LOCATE_INFO_IN_INDEX_HEADERID, NULL)))
    return 0;
  if (fread(headerData, sizeof (*headerData),
            1, fp) != 1)
    return 0;
  return 1;
}

struct seqRevMapEntry
{
  Seqpos bwtPos, origPos;
};

struct addLocateInfoState
{
  Seqpos seqLen;
  void *spReadState, *origSeqState;
  const MRAEnc *alphabet;
  int *specialRanges;
  SeqposReadFunc readSeqpos;
  GetOrigSeqSym readOrigSeq;
  unsigned locateInterval, bitsPerOrigPos, bitsPerSeqpos;
  int featureToggles;
  size_t revMapCacheSize;
  struct seqRevMapEntry *revMapCache;
};

static inline unsigned
bitsPerPosUpperBound(struct addLocateInfoState *state)
{
  unsigned bitsPerPos = 0;
  if (state->featureToggles & BWTLocateCount)
    /* +1 to account for the contribution to count value (which will
     * actually use log(blockLen) instead of blockLen bits, i.e. less
     * than 1 bit per position) */
    bitsPerPos += state->bitsPerSeqpos + 1;
  bitsPerPos += state->bitsPerOrigPos; /* stored for every position */
  return bitsPerPos;
}

static void
initAddLocateInfoState(struct addLocateInfoState *state,
                       GetOrigSeqSym readOrigSeq, void *origSeqState,
                       SeqposReadFunc readSeqpos, void *spReadState,
                       const MRAEnc *alphabet, int *specialRanges,
                       Seqpos srcLen, const struct bwtParam *params)
{
  Seqpos lastPos;
  unsigned aggregationExpVal;
  unsigned locateInterval;
  assert(state);
  state->alphabet = alphabet;
  state->seqLen = srcLen;
  state->specialRanges = specialRanges;
  state->readSeqpos = readSeqpos;
  state->spReadState = spReadState;
  state->readOrigSeq = readOrigSeq;
  state->origSeqState = origSeqState;
  state->featureToggles = params->featureToggles;
  aggregationExpVal = estimateSegmentSize(&params->seqParams,
                                          params->baseType);
  locateInterval = params->locateInterval;
  lastPos = srcLen - 1;
  state->locateInterval = locateInterval;
  state->bitsPerSeqpos = requiredSeqposBits(lastPos);
  if (locateInterval)
  {
    if (params->featureToggles & BWTProperlySorted)
      state->bitsPerOrigPos = requiredSeqposBits(lastPos/locateInterval);
    else
      state->bitsPerOrigPos = requiredSeqposBits(lastPos);
    state->revMapCacheSize = aggregationExpVal;
    state->revMapCache = ma_malloc(sizeof (state->revMapCache[0])
                                   * state->revMapCacheSize);
  }
  else
  {
    state->bitsPerOrigPos = 0;
    state->revMapCacheSize = 0;
    state->revMapCache = NULL;
  }
}

static void
destructAddLocateInfoState(struct addLocateInfoState *state)
{
  ma_free(state->revMapCache);
}

static BitOffset
addLocateInfo(BitString cwDest, BitOffset cwOffset,
              BitString varDest, BitOffset varOffset,
              Seqpos start, Seqpos len, void *cbState, Error *err)
{
  BitOffset bitsWritten = 0;
  struct addLocateInfoState *state = cbState;
  unsigned bitsPerBWTPos, bitsPerOrigPos;
  int retcode;
  assert(varDest && cbState);
  /* 0. resize caches if necessary */
  if (len > state->revMapCacheSize)
  {
    state->revMapCache = ma_realloc(state->revMapCache,
                                    len * sizeof (state->revMapCache[0]));
    state->revMapCacheSize = len;
  }
  bitsPerBWTPos = requiredSeqposBits(len - 1);
  bitsPerOrigPos = state->bitsPerOrigPos;
  /* 1. read rangeLen suffix array indices from suftab */
  {
    Seqpos i, mapVal;
    size_t revMapCacheLen = 0;
    int properlySorted = state->featureToggles & BWTProperlySorted,
      locateBitmap = state->featureToggles & BWTLocateBitmap;
    for (i = 0; i < len; ++i)
    {
      int specialRegularSymTransition = 0;
      /* 1.a read array index*/
      if ((retcode = state->readSeqpos(state->spReadState, &mapVal, 1, err))
          != 1)
        return (BitOffset)-1;
      /* 1.b find current symbol and compare to special ranges */
      if (!properlySorted)
      {
        Symbol syms[2];
        if (mapVal > 0 && mapVal < state->seqLen - 1)
        {
          if ((retcode = state->readOrigSeq(
                 state->origSeqState, syms, mapVal - 1, 2)) != 2)
            return (BitOffset)-1;
        }
        else if (mapVal == 0)
        {
          syms[0] = UNDEFBWTCHAR;
          if ((retcode = state->readOrigSeq(
                 state->origSeqState, syms, 0, 1))!= 1)
            return (BitOffset)-1;
        }
        else /* mapVal == seqLen - 1 */
        {
          if ((retcode = state->readOrigSeq(state->origSeqState,
                                            syms, state->seqLen - 2, 1)) != 1)
            return (BitOffset)-1;
          syms[1] = UNDEFBWTCHAR;
        }
        if (MRAEncSymbolIsInSelectedRanges(
              state->alphabet, MRAEncMapSymbol(state->alphabet, syms[0]),
              SPECIAL_RANGE, state->specialRanges)
            && MRAEncSymbolIsInSelectedRanges(
              state->alphabet, MRAEncMapSymbol(state->alphabet, syms[1]),
              NORMAL_RANGE, state->specialRanges))
          specialRegularSymTransition = 1;
        else
          specialRegularSymTransition = 0;
      }
      /* 1.c check wether the index into the original sequence is an
       * even multiple of the sampling interval, or we have just
       * walked from a special range to a regular symbol */
      if (!(mapVal % state->locateInterval)
          || (!properlySorted && specialRegularSymTransition))
      {
        /* 1.c.1 enter index into cache */
        state->revMapCache[revMapCacheLen].bwtPos = i;
        state->revMapCache[revMapCacheLen].origPos
          = properlySorted?mapVal / state->locateInterval:mapVal;
        ++revMapCacheLen;
        /* 1.c.2 mark position in bwt sequence */
        if (locateBitmap)
          bsSetBit(cwDest, cwOffset + i);
      }
      else
        if (locateBitmap)
          bsClearBit(cwDest, cwOffset + i);
    }
    /* 2. copy revMapCache into output */
    {
      int locateCount = state->featureToggles & BWTLocateCount;
      if (locateCount)
      {
        unsigned bitsPerCount = requiredSeqposBits(len);
        bsStoreSeqpos(varDest, varOffset + bitsWritten, bitsPerCount,
                      revMapCacheLen);
        bitsWritten += bitsPerCount;
      }
      for (i = 0; i < revMapCacheLen; ++i)
      {
        if (locateCount)
        {
          bsStoreSeqpos(varDest, varOffset + bitsWritten, bitsPerBWTPos,
                        state->revMapCache[i].bwtPos);
          bitsWritten += bitsPerBWTPos;
        }
        bsStoreSeqpos(varDest, varOffset + bitsWritten, bitsPerOrigPos,
                      state->revMapCache[i].origPos);
        bitsWritten += bitsPerOrigPos;
      }
    }
  }
  return bitsWritten;
}

extern EISeq *
createBWTSeqGeneric(const struct bwtParam *params,
                    indexCreateFunc createIndex, void *baseSrc,
                    Seqpos totalLen,
                    const MRAEnc *alphabet, int *specialRanges,
                    GetOrigSeqSym readOrigSeq, void *origSeqState,
                    SeqposReadFunc readNextSeqpos, void *spReadState,
                    reportLongest lrepFunc, void *lrepState, Error *err)
{
  struct encIdxSeq *baseSeqIdx = NULL;
  struct addLocateInfoState varState;
  bool varStateIsInitialized = false;
  unsigned locateInterval;
  assert(baseSrc && params && err);
  error_check(err);
  locateInterval = params->locateInterval;
  do
  {
    struct locateHeaderWriteInfo headerData
      = { locateInterval, lrepFunc, lrepState, params->featureToggles };
    void *p[] = { &headerData };
    uint16_t headerIDs[] = { LOCATE_INFO_IN_INDEX_HEADERID };
    uint32_t headerSizes[] = { LOCATE_HEADER_SIZE };
    headerWriteFunc headerFuncs[] = { writeLocateInfoHeader };
    error_unset(err);
    initAddLocateInfoState(&varState,
                           readOrigSeq, origSeqState,
                           readNextSeqpos, spReadState,
                           alphabet, specialRanges,
                           totalLen, params);
    varStateIsInitialized = true;
    if (locateInterval)
    {
      if (!(baseSeqIdx
            = createIndex(baseSrc, totalLen, params->projectName,
                          &params->seqParams, sizeof (p)/sizeof (p[0]),
                          headerIDs, headerSizes, headerFuncs, p,
                          addLocateInfo,
                          /* one bit per position if using bitmap */
                          (params->featureToggles & BWTLocateBitmap)?1:0,
                          bitsPerPosUpperBound(&varState),
                          &varState, err)))
        break;
    }
    else
    {
      if (!(baseSeqIdx
            = createIndex(baseSrc, totalLen, params->projectName,
                          &params->seqParams, 0, NULL, NULL, NULL,
                          NULL, NULL, 0, 0, &varState, err)))
        break;
    }
  } while (0);
  if (varStateIsInitialized)
    destructAddLocateInfoState(&varState);
  return baseSeqIdx;
}

static inline BitOffset
searchLocateCountMark(const BWTSeq *bwtSeq, Seqpos pos,
                      struct extBitsRetrieval *extBits)
{
  unsigned i, numMarks, bitsPerCount;
  BitOffset markOffset;
  EISRetrieveExtraBits(bwtSeq->seqIdx, pos, EBRF_RETRIEVE_CWBITS
                       | EBRF_RETRIEVE_VARBITS, extBits, bwtSeq->hint);
  markOffset = extBits->varOffset;
  bitsPerCount = requiredSeqposBits(extBits->len);
  numMarks = bsGetSeqpos(extBits->varPart, markOffset, bitsPerCount);
  if (numMarks)
  {
    unsigned bitsPerBWTPos = requiredSeqposBits(extBits->len - 1),
      bitsPerOrigPos = requiredSeqposBits(BWTSeqLength(bwtSeq) - 1);
    Seqpos cmpPos = pos - extBits->start;
    markOffset += bitsPerCount;
    for (i = 0; i < numMarks; ++i)
    {
      Seqpos markedPos = bsGetSeqpos(extBits->varPart, markOffset,
                                     bitsPerBWTPos);
      if (markedPos < cmpPos)
        markOffset += bitsPerBWTPos + bitsPerOrigPos;
      else if (markedPos > cmpPos)
        break;
      else if (markedPos == cmpPos)
        return markOffset + bitsPerBWTPos;
    }
  }
  return 0;
}

extern int
BWTSeqPosHasLocateInfo(const BWTSeq *bwtSeq, Seqpos pos,
                       struct extBitsRetrieval *extBits)
{
  if (bwtSeq->featureToggles & BWTLocateBitmap)
  {
    EISRetrieveExtraBits(bwtSeq->seqIdx, pos, EBRF_RETRIEVE_CWBITS, extBits,
                         bwtSeq->hint);
    return bsGetBit(extBits->cwPart, extBits->cwOffset + pos - extBits->start);
  }
  else if (bwtSeq->featureToggles & BWTLocateCount)
  {
    BitOffset markOffset = searchLocateCountMark(bwtSeq, pos, extBits);
    return markOffset != 0;
  }
#ifndef NDEBUG
  else
  {
    fputs("Trying to locate in BWT sequence index without locate information.",
          stderr);
    abort();
  }
#endif
  return 0;
}

extern Seqpos
BWTSeqLocateMatch(const BWTSeq *bwtSeq, Seqpos pos,
                  struct extBitsRetrieval *extBits)
{
  if (bwtSeq->featureToggles & BWTLocateBitmap)
  {
    Seqpos nextLocate = pos;
    unsigned locateOffset = 0;
    while (!BWTSeqPosHasLocateInfo(bwtSeq, nextLocate, extBits))
      nextLocate = BWTSeqLFMap(bwtSeq, nextLocate), ++locateOffset;
    EISRetrieveExtraBits(bwtSeq->seqIdx, nextLocate,
                         EBRF_RETRIEVE_CWBITS | EBRF_RETRIEVE_VARBITS,
                         extBits, bwtSeq->hint);
    {
      Seqpos maxPosVal = ((bwtSeq->featureToggles & BWTProperlySorted)?
                          (BWTSeqLength(bwtSeq) - 1)
                          /bwtSeq->locateSampleInterval:
                          BWTSeqLength(bwtSeq) - 1);
      unsigned bitsPerCount = requiredSeqposBits(extBits->len),
        bitsPerBWTPos = requiredSeqposBits(extBits->len - 1),
        bitsPerOrigPos = requiredSeqposBits(maxPosVal);
      BitOffset locateRecordIndex =
        bs1BitsCount(extBits->cwPart, extBits->cwOffset,
                     nextLocate - extBits->start),
        locateRecordOffset = ((bwtSeq->featureToggles & BWTLocateCount?
                               bitsPerBWTPos:0) + bitsPerOrigPos)
        * locateRecordIndex
        + ((bwtSeq->featureToggles & BWTLocateCount)?bitsPerCount:0);
      Seqpos matchPos =
        bsGetSeqpos(
          extBits->varPart, extBits->varOffset + locateRecordOffset
          + ((bwtSeq->featureToggles & BWTLocateCount)?bitsPerBWTPos:0),
          bitsPerOrigPos);
      if (bwtSeq->featureToggles & BWTProperlySorted)
        matchPos = matchPos * bwtSeq->locateSampleInterval;
      matchPos += locateOffset;
      assert(!(bwtSeq->featureToggles & BWTLocateCount)
             || bsGetSeqpos(extBits->varPart,
                            extBits->varOffset + locateRecordOffset,
                            bitsPerBWTPos)
             == nextLocate - extBits->start);
      return matchPos;
    }
  }
  else if (bwtSeq->featureToggles & BWTLocateCount)
  {
    BitOffset markOffset;
    Seqpos nextLocate = pos, matchPos;
    unsigned locateOffset = 0;  /* mark is at most locateInterval
                                 * positions away */
    unsigned bitsPerOrigPos
      = requiredSeqposBits(BWTSeqLength(bwtSeq) - 1);
    while ((markOffset = searchLocateCountMark(bwtSeq, nextLocate,
                                               extBits)) == 0)
      nextLocate = BWTSeqLFMap(bwtSeq, nextLocate), ++locateOffset;
    matchPos = bsGetSeqpos(extBits->varPart, markOffset, bitsPerOrigPos)
      + locateOffset;
    return matchPos;
  }
  /* Internal error: Trying to locate in BWT sequence index without locate
     information. */
   assert(0);
   return 0; /* shut up compiler */
}
