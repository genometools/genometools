/*
  Copyright (C) 2007,2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#include "core/chardef.h"
#include "core/log.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "match/eis-specialsrank.h"
#include "match/eis-bitpackseqpos.h"
#include "match/eis-bwtseq.h"
#include "match/eis-bwtseq-extinfo.h"
#include "match/eis-bwtseq-priv.h"
#include "match/eis-bwtseq-context.h"
#include "match/eis-headerid.h"
#include "match/eis-mrangealphabet.h"
#include "match/eis-sa-common.h"

/**
 * @file eis-bwtseq-extinfo.c generic methods for bwt index creation and
 * everything related to storing/retrieving locate information
 */

struct locateHeader
{
  unsigned long rot0Pos;
  unsigned locateInterval;
  int featureToggles;
};

enum {
  LOCATE_HEADER_SIZE = sizeof (struct locateHeader),
};

struct locateHeaderWriteInfo
{
  SASeqSrc *src;
  unsigned locateInterval;
  int featureToggles;
};

static int
writeLocateInfoHeader(FILE *fp, void *cbData)
{
  struct locateHeader headerData;
  Definedunsignedlong rot0Pos;
  const struct locateHeaderWriteInfo *headerSrc = cbData;
  gt_assert(cbData);
  headerData.locateInterval = headerSrc->locateInterval;
  rot0Pos = SASSGetRot0Pos(headerSrc->src);
  if (!rot0Pos.defined)
  {
    fputs("Invalid index construction: position of suffix 0 unknown!\n",
          stderr);
    abort();
  }
  headerData.rot0Pos = rot0Pos.valueunsignedlong;
  headerData.featureToggles = headerSrc->featureToggles;
  gt_xfwrite(&headerData, sizeof (headerData), 1, fp);
  return 1;
}

/**
 * @return 0 on error
 */
static inline int
readLocateInfoHeader(EISeq *seqIdx, struct locateHeader *headerData)
{
  FILE *fp;
  gt_assert(seqIdx && headerData);
  if (!(fp = EISSeekToHeader(seqIdx, LOCATE_INFO_IN_INDEX_HEADERID, NULL)))
    return 0;
  if (fread(headerData, sizeof (*headerData), 1, fp) != 1)
    return 0;
  return 1;
}

struct sortModeHeader
{
  uint32_t bitsPerOrigRank;
  const MRAEnc *alphabet;
  const enum rangeSortMode *rangeSort;
};

static inline uint32_t
computeSortModeHeaderSize(const MRAEnc *alphabet)
{
  return sizeof (uint32_t) + sizeof (int16_t) * MRAEncGetNumRanges(alphabet);
}

static int
writeRankSortHeader(FILE *fp, void *cbData)
{
  struct sortModeHeader *headerData = cbData;
  gt_assert(cbData);
  gt_xfwrite(&headerData->bitsPerOrigRank, sizeof (headerData->bitsPerOrigRank),
             1, fp);
  {
    size_t i, numRanges = MRAEncGetNumRanges(headerData->alphabet);
    for (i = 0; i < numRanges; ++i)
    {
      int16_t mode = headerData->rangeSort[i];
      gt_xfwrite(&mode, sizeof (mode), 1, fp);
    }
  }
  return 1;
}

static inline int
readRankSortHeader(EISeq *seqIdx, uint32_t *bitsPerOrigRank,
                   const MRAEnc *alphabet,
                   enum rangeSortMode *rangeSort)
{
  FILE *fp;
  gt_assert(seqIdx && alphabet && bitsPerOrigRank && rangeSort);
  if (!(fp = EISSeekToHeader(seqIdx, RANK_SORT_HEADERID, NULL)))
    return 0;
  if (fread(bitsPerOrigRank, sizeof (*bitsPerOrigRank), 1, fp) != 1)
    return 0;
  {
    uint16_t mode;
    size_t i, numRanges = MRAEncGetNumRanges(alphabet);
    for (i = 0; i < numRanges; ++i)
    {
      if (fread(&mode, sizeof (mode), 1, fp) != 1)
        return 0;
      rangeSort[i] = mode;
    }
  }
  return 1;
}

struct seqRevMapEntry
{
  unsigned long bwtPos, origPos;
};

struct addLocateInfoState
{
  unsigned long seqLen, extraLocMarksUpperBound;
  const MRAEnc *alphabet;
  const enum rangeSortMode *rangeSort;
  SeqDataReader readUlong;
  RandomSeqAccessor origSeqAccess;
  unsigned locateInterval, bitsPerOrigPos, bitsPerUlong,
    bitsPerOrigRank;
  const SpecialsRankLookup *sprTable;
  int featureToggles;
  size_t revMapQueueSize;
  struct seqRevMapEntry *revMapQueue;
  size_t origRanksQueueSize;
  unsigned long *origRanksQueue;
  BWTSeqContextRetrieverFactory *ctxFactory;
};

static inline unsigned
bitsPerPosUpperBoundWithoutSegmentAggregates(struct addLocateInfoState *state)
{
  unsigned bitsPerPos = 0;
  if (state->featureToggles & BWTLocateCount)
    /* bitsPerUlong is actually an upper bound of log(segmentLen) to
     * identify each locate mark */
    bitsPerPos += state->bitsPerUlong;
  bitsPerPos += state->bitsPerOrigPos; /* stored for every position */
  /* rank sorting data if needed */
  bitsPerPos += state->bitsPerOrigRank;
  return bitsPerPos;
}

static inline unsigned
bitsPerPosUpperBound(struct addLocateInfoState *state)
{
  unsigned bitsPerPos = bitsPerPosUpperBoundWithoutSegmentAggregates(state);
  if (state->featureToggles & BWTLocateCount)
    /* +1 to account for the contribution to count value (which will
     * actually use log(blockLen) instead of blockLen bits, i.e. less
     * than 1 bit per position) */
    bitsPerPos += 1;
  return bitsPerPos;
}

static bool
locBitsUpperBounds(void *cbState, struct segmentDesc *desc,
                   size_t numSegmentDesc, struct varBitsEstimate *result)
{
  struct addLocateInfoState *state = cbState;
  gt_assert(cbState);
  if (state->featureToggles & (BWTLocateCount | BWTLocateBitmap))
  {
    unsigned maxBitsPerPos = bitsPerPosUpperBound(state);
    result->maxBitsPerPos = maxBitsPerPos;
    if (desc)
    {
      size_t i, maxSegLen = 0;
      BitOffset maxBitsTotal = 0;
      unsigned long numSegmentsTotal = 0;
      unsigned bitsPerOrigRank = state->bitsPerOrigRank;
      for (i = 0; i < numSegmentDesc; ++i)
      {
        size_t len = desc[i].len;
        maxSegLen = MAX(len, maxSegLen);
        if (state->featureToggles & BWTLocateCount)
          maxBitsTotal += requiredUlongBits(len) * desc[i].repeatCount;
        numSegmentsTotal += desc[i].repeatCount;
      }
      maxBitsTotal +=
        (state->seqLen / state->locateInterval + state->extraLocMarksUpperBound)
        * (((state->featureToggles & BWTLocateCount) ?
            requiredUlongBits(maxSegLen) : 0)
           + state->bitsPerOrigPos);
      if (bitsPerOrigRank)
      {
        unsigned long maxRank = specialsRank(state->sprTable, state->seqLen);
        /* rank values stored */
        maxBitsTotal += maxRank * bitsPerOrigRank;
      }
      result->maxBitsPerBucket =
        /* locate marks and rank sorts */
        maxSegLen * bitsPerPosUpperBoundWithoutSegmentAggregates(state)
        /* count of locate marks */
        + ((state->featureToggles & BWTLocateCount) ?
           requiredUlongBits(maxSegLen) : 0);
      result->maxBitsTotal = maxBitsTotal;
      return true;
    }
    else
      return false;
  }
  else
  {
    result->maxBitsTotal = result->maxBitsPerBucket
      = result->maxBitsPerPos = 0;
    return false;
  }
}

static void
initAddLocateInfoState(struct addLocateInfoState *state,
                       RandomSeqAccessor origSeqAccess,
                       SeqDataReader readUlong,
                       const MRAEnc *alphabet, const struct seqStats *stats,
                       const enum rangeSortMode *rangeSort,
                       unsigned long srcLen, const struct bwtParam *params,
                       const SpecialsRankLookup *sprTable,
                       unsigned bitsPerOrigRank,
                       BWTSeqContextRetrieverFactory *ctxFactory)
{
  unsigned long lastPos;
  unsigned aggregationExpVal;
  unsigned locateInterval;
  gt_assert(state);
  state->alphabet = alphabet;
  state->seqLen = srcLen;
  state->rangeSort = rangeSort;
  state->readUlong = readUlong;
  state->origSeqAccess = origSeqAccess;
  state->featureToggles = params->featureToggles;
  aggregationExpVal = gt_estimateSegmentSize(&params->seqParams);
  locateInterval = params->locateInterval;
  lastPos = srcLen - 1;
  state->locateInterval = locateInterval;
  state->bitsPerUlong = requiredUlongBits(lastPos);
  state->bitsPerOrigRank = bitsPerOrigRank;
  state->sprTable = sprTable;
  if (locateInterval)
  {
    if (params->featureToggles & BWTReversiblySorted)
      state->bitsPerOrigPos = requiredUlongBits(lastPos/locateInterval);
    else
      state->bitsPerOrigPos = requiredUlongBits(lastPos);
    state->origRanksQueueSize = state->revMapQueueSize =aggregationExpVal;
    state->revMapQueue = gt_malloc(sizeof (state->revMapQueue[0])
                                   * state->revMapQueueSize);
    state->origRanksQueue = gt_malloc(sizeof (state->origRanksQueue[0])
                                      * state->origRanksQueueSize);
    state->bitsPerOrigRank = bitsPerOrigRank;
    /* if there are rank sorted symbols, but no rank sorting is
     * available, extra locate marks will be inserted, at most every
     * second position can be marked this way because for every
     * SORTMODE_VALUE -> SORTMODE_{RANk|UNDEF} transition there must be one
     * inverse transition */
    if (!(params->featureToggles & BWTReversiblySorted)
        && locateInterval > 1)
    {
      unsigned long stdLocMarks = srcLen / locateInterval,
        extraLocMarksUpperBound = MIN(srcLen/2, srcLen - stdLocMarks);
      if (stats)
      {
        unsigned long nonValSortSyms = 0;
        unsigned i;
        for (i = 0; i <= UINT8_MAX; ++i)
          if (MRAEncSymbolHasValidMapping(alphabet, i)
              && !gt_MRAEncSymbolIsInSelectedRanges(
                alphabet, MRAEncMapSymbol(alphabet, i),
                SORTMODE_VALUE, (int *)rangeSort))
            nonValSortSyms += stats->symbolDistributionTable[i];
        extraLocMarksUpperBound = MIN3(extraLocMarksUpperBound, nonValSortSyms,
                                       srcLen - nonValSortSyms);
      }
      state->extraLocMarksUpperBound = extraLocMarksUpperBound;
    }
    else
    {
      state->extraLocMarksUpperBound = 0;
    }
  }
  else
  {
    state->extraLocMarksUpperBound = 0;
    state->bitsPerOrigPos = state->bitsPerOrigRank = 0;
    state->origRanksQueueSize = state->revMapQueueSize = 0;
    state->origRanksQueue = NULL;
    state->revMapQueue = NULL;
  }
  state->ctxFactory = ctxFactory;
}

static void
destructAddLocateInfoState(struct addLocateInfoState *state)
{
  gt_free(state->revMapQueue);
  gt_free(state->origRanksQueue);
}

static inline int
isSortModeTransition(RandomSeqAccessor origSeqAccess, unsigned long seqLen,
                     const MRAEnc *alphabet,
                     const enum rangeSortMode *rangeSort,
                     unsigned long pos)
{
  Symbol syms[2];
  if (pos > 0 && pos < seqLen - 1)
  {
#ifndef NDEBUG
    int retcode =
#endif
      accessSequence(origSeqAccess, syms, pos - 1, 2);
    gt_assert(retcode == 2);
  }
  else if (pos == 0)
  {
#ifndef NDEBUG
    int retcode =
#endif
      accessSequence(origSeqAccess, syms + 1, 0, 1);
    gt_assert(retcode == 1);
    syms[0] = UNDEFBWTCHAR;
  }
  else /* pos == seqLen - 1 */
  {
#ifndef NDEBUG
    int retcode =
#endif
      accessSequence(origSeqAccess, syms, seqLen - 2, 1);
    gt_assert(retcode == 1);
    syms[1] = UNDEFBWTCHAR;
  }
  {
    AlphabetRangeID range[2];
    range[0] = MRAEncGetRangeOfSymbol(alphabet,
                                      MRAEncMapSymbol(alphabet, syms[0]));
    range[1] = MRAEncGetRangeOfSymbol(alphabet,
                                      MRAEncMapSymbol(alphabet, syms[1]));
    return rangeSort[range[0]] != rangeSort[range[1]];
  }
}

static BitOffset
addLocateInfo(BitString cwDest, BitOffset cwOffset,
              BitString varDest, BitOffset varOffset,
              GT_UNUSED unsigned long start, unsigned long len, void *cbState)
{
  BitOffset bitsWritten = 0;
  struct addLocateInfoState *state = cbState;
  unsigned bitsPerBWTPos, bitsPerOrigPos, bitsPerOrigRank;
  int retcode;
  unsigned locateInterval;
  gt_assert(varDest && cbState);
  locateInterval = state->locateInterval;
  /* 0. resize caches if necessary */
  if (locateInterval && len > state->revMapQueueSize)
  {
    state->revMapQueue = gt_realloc(state->revMapQueue,
                                    len * sizeof (state->revMapQueue[0]));
    state->revMapQueueSize = len;
  }
  if (state->sprTable && len > state->origRanksQueueSize)
  {
    state->origRanksQueue = gt_realloc(state->origRanksQueue,
                                          len *
                                          sizeof (state->origRanksQueue[0]));
    state->origRanksQueueSize = len;
  }
  bitsPerBWTPos = requiredUlongBits(len - 1);
  bitsPerOrigPos = state->bitsPerOrigPos;
  bitsPerOrigRank = state->bitsPerOrigRank;
  {
    unsigned long i, mapVal;
    size_t revMapQueueLen = 0, origRanksQueueLen = 0;
    int reversiblySorted = state->featureToggles & BWTReversiblySorted,
      locateBitmap = state->featureToggles & BWTLocateBitmap;
    /* read len suffix array indices from suftab */
    /* TODO: decide  on reasonable buffering of mapVal */
    for (i = 0; i < len; ++i)
    {
      /* add locate data if necessary */
      if (locateInterval)
      {
        int insertExtraLocateMark = 0;
        /* 1.a read array index*/
        if ((retcode = SDRRead(state->readUlong, &mapVal, 1)) != 1)
          return (BitOffset)-1;
        /* 1.b find current symbol and compare to special ranges */
        insertExtraLocateMark = !reversiblySorted &&
          isSortModeTransition(state->origSeqAccess, state->seqLen,
                               state->alphabet, state->rangeSort,
                               mapVal);
        /* 1.c check wether the index into the original sequence is an
         * even multiple of the sampling interval, or a not-reversible
         * sort mode transition occurred */
        if (!(mapVal % locateInterval) || insertExtraLocateMark)
        {
          /* 1.c.1 enter index into cache */
          state->revMapQueue[revMapQueueLen].bwtPos = i;
          state->revMapQueue[revMapQueueLen].origPos
            = reversiblySorted ? mapVal/state->locateInterval : mapVal;
          ++revMapQueueLen;
          /* 1.c.2 mark position in bwt sequence */
          if (locateBitmap)
            bsSetBit(cwDest, cwOffset + i);
        }
        else
          if (locateBitmap)
            gt_bsClearBit(cwDest, cwOffset + i);
      }
      /* and add data for extra sort information */
      if (bitsPerOrigRank)
      {
        Symbol BWTSym;
        const MRAEnc *alphabet = state->alphabet;
        if (mapVal != 0)
          accessSequence(state->origSeqAccess, &BWTSym, mapVal - 1, 1);
        else
          BWTSym = UNDEFBWTCHAR;
        if (state->rangeSort[MRAEncGetRangeOfSymbol(
                alphabet, MRAEncMapSymbol(alphabet, BWTSym))]
            == SORTMODE_RANK)
        {
          state->origRanksQueue[origRanksQueueLen++] =
            specialsRank(state->sprTable,
                         mapVal != 0 ? mapVal - 1 : state->seqLen - 1);
        }
      }
      if (state->ctxFactory)
      {
        gt_BWTSCRFMapAdvance(state->ctxFactory, &mapVal, 1);
      }
    }
    /* 2. copy revMapQueue into output */
    if (locateInterval)
    {
      int locateCount = state->featureToggles & BWTLocateCount;
      if (locateCount)
      {
        unsigned bitsPerCount = requiredUlongBits(len);
        gt_bsStoreUlong(varDest, varOffset + bitsWritten, bitsPerCount,
                      revMapQueueLen);
        bitsWritten += bitsPerCount;
      }
      for (i = 0; i < revMapQueueLen; ++i)
      {
        if (locateCount)
        {
          gt_bsStoreUlong(varDest, varOffset + bitsWritten, bitsPerBWTPos,
                        state->revMapQueue[i].bwtPos);
          bitsWritten += bitsPerBWTPos;
        }
        gt_bsStoreUlong(varDest, varOffset + bitsWritten, bitsPerOrigPos,
                      state->revMapQueue[i].origPos);
        bitsWritten += bitsPerOrigPos;
      }
    }
    if (bitsPerOrigRank)
    {
      unsigned bitsPerOrigRank = state->bitsPerOrigRank;
      if (origRanksQueueLen)
      {
        gt_bsStoreUniformUlongArray(varDest, varOffset + bitsWritten,
                                  bitsPerOrigRank, origRanksQueueLen,
#ifdef _LP64
                             (uint64_t*) state->origRanksQueue);
#else
                             (uint32_t*) state->origRanksQueue);
#endif
        bitsWritten += bitsPerOrigRank * origRanksQueueLen;
      }
    }
  }
  return bitsWritten;
}

static inline int
sortModeHeaderNeeded(const MRAEnc *alphabet,
                     const enum rangeSortMode *rangeSort,
                     const SpecialsRankLookup *sprTable)
{
  bool hasRankSortedRanges = false;
  AlphabetRangeID i, numRanges = MRAEncGetNumRanges(alphabet);
  for (i = 0; i < numRanges; ++i)
    hasRankSortedRanges |= (rangeSort[i] == SORTMODE_RANK?1:0);
  return (hasRankSortedRanges && sprTable);
}

EISeq *
gt_createBWTSeqGeneric(const struct bwtParam *params,
                       indexCreateFunc createIndex,
                       SASeqSrc *src,
                       const enum rangeSortMode rangeSort[],
                       const SpecialsRankLookup *sprTable,
                       GtError *err)
{
  struct encIdxSeq *baseSeqIdx = NULL;
  struct addLocateInfoState varState;
  bool varStateIsInitialized = false;
  unsigned locateInterval;
  BWTSeqContextRetrieverFactory *buildContextMap = NULL;
  gt_assert(src && params && err);
  gt_error_check(err);
  locateInterval = params->locateInterval;
  do
  {
    struct locateHeaderWriteInfo locHeaderData
      = { src, locateInterval, params->featureToggles };
    struct sortModeHeader sortModeHeader;
    void *p[] = { &locHeaderData , &sortModeHeader };
    uint16_t headerIDs[] = { LOCATE_INFO_IN_INDEX_HEADERID,
                             RANK_SORT_HEADERID };
    uint32_t headerSizes[] = { LOCATE_HEADER_SIZE,
                               0 };
    headerWriteFunc headerFuncs[] = { writeLocateInfoHeader,
                                      writeRankSortHeader };
    size_t numHeaders = 0;
    unsigned bitsPerOrigRank = 0;
    unsigned long totalLen = SASSGetLength(src);
    const MRAEnc *alphabet = SASSGetMRAEnc(src);
    MRAEnc *baseAlphabet = SASSNewMRAEnc(src);
    /* FIXME: this  has to work also when locateInterval == 0 and
     * sprTable != NULL */
    if (params->ctxMapILog != CTX_MAP_ILOG_NOMAP)
      buildContextMap = gt_newBWTSeqContextRetrieverFactory(totalLen,
                                                         params->ctxMapILog);
    if (locateInterval)
    {
      ++numHeaders;
      if (sortModeHeaderNeeded(alphabet, rangeSort, sprTable))
      {
        unsigned long
#ifndef NDEBUG
          origSeqLen = gt_encseq_total_length(
                                                gt_SPRTGetOrigEncseq(sprTable)),
#endif
          maxRank;
        gt_assert(origSeqLen == totalLen - 1);
        maxRank = specialsRank(sprTable, totalLen - 1);
        bitsPerOrigRank = sortModeHeader.bitsPerOrigRank
          = requiredUlongBits(maxRank);
        sortModeHeader.alphabet = alphabet;
        sortModeHeader.rangeSort = rangeSort;
        headerSizes[1] = computeSortModeHeaderSize(alphabet);
        ++numHeaders;
      }
      {
        SeqDataReader readSfxIdx = SASSCreateReader(src, SFX_REQUEST_SUFTAB);
        if (SDRIsValid(readSfxIdx))
        {
          initAddLocateInfoState(
            &varState, SASSGetOrigSeqAccessor(src), readSfxIdx,
            alphabet, SASSGetSeqStats(src), rangeSort, totalLen, params,
            bitsPerOrigRank?sprTable:NULL, bitsPerOrigRank,
            buildContextMap);
          varStateIsInitialized = true;
        }
        else
        {
          gt_error_set(err, "error: locate sampling requested but not available"
                       " for project %s\n", gt_str_get(params->projectName));
        }
      }
    }
    if (!(baseSeqIdx
          = createIndex(totalLen, gt_str_get(params->projectName), baseAlphabet,
                        SASSGetSeqStats(src),
                        SASSCreateReader(src, SFX_REQUEST_BWTTAB),
                        &params->seqParams, numHeaders,
                        headerIDs, headerSizes, headerFuncs, p,
                        locateInterval?addLocateInfo:NULL,
                        /* one bit per position if using bitmap */
                        (params->featureToggles & BWTLocateBitmap)?1:0,
                        locateInterval?locBitsUpperBounds:NULL, &varState,
                        err)))
      break;
    if (buildContextMap)
    {
      if (!gt_BWTSCRFFinished(buildContextMap))
      {
        fputs("error: context table construction incomplete!\n", stderr);
      }
      else
      {
        BWTSeqContextRetriever *ctxRetrieve =
          gt_BWTSCRFGet(buildContextMap, NULL, gt_str_get(params->projectName));
        gt_deleteBWTSeqCR(ctxRetrieve);
      }
    }
  } while (0);
  if (buildContextMap)
    gt_deleteBWTSeqContextRetrieverFactory(buildContextMap);
  if (varStateIsInitialized)
    destructAddLocateInfoState(&varState);
  return baseSeqIdx;
}

static inline BitOffset
searchLocateCountMark(const BWTSeq *bwtSeq, unsigned long pos,
                      struct extBitsRetrieval *extBits)
{
  unsigned i, numMarks, bitsPerCount;
  BitOffset markOffset;
  EISRetrieveExtraBits(bwtSeq->seqIdx, pos, EBRF_RETRIEVE_CWBITS
                       | EBRF_RETRIEVE_VARBITS, extBits, bwtSeq->hint);
  markOffset = extBits->varOffset;
  bitsPerCount = requiredUlongBits(extBits->len);
  numMarks = gt_bsGetUlong(extBits->varPart, markOffset, bitsPerCount);
  if (numMarks)
  {
    unsigned bitsPerBWTPos = requiredUlongBits(extBits->len - 1),
      bitsPerOrigPos = requiredUlongBits(
        ((bwtSeq->featureToggles & BWTReversiblySorted) ?
         (BWTSeqLength(bwtSeq) - 1) / bwtSeq->locateSampleInterval :
         BWTSeqLength(bwtSeq) - 1));
    unsigned long cmpPos = pos - extBits->start;
    markOffset += bitsPerCount;
    for (i = 0; i < numMarks; ++i)
    {
      unsigned long markedPos = gt_bsGetUlong(extBits->varPart, markOffset,
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

int
gt_BWTSeqPosHasLocateInfo(const BWTSeq *bwtSeq, unsigned long pos,
                       struct extBitsRetrieval *extBits)
{
  if (bwtSeq->featureToggles & BWTLocateBitmap)
  {
    EISRetrieveExtraBits(bwtSeq->seqIdx, pos, EBRF_RETRIEVE_CWBITS, extBits,
                         bwtSeq->hint);
    return gt_bsGetBit(extBits->cwPart,
                       extBits->cwOffset + pos - extBits->start);
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

unsigned long
gt_BWTSeqLocateMatch(const BWTSeq *bwtSeq, unsigned long pos,
                  struct extBitsRetrieval *extBits)
{
  if (bwtSeq->featureToggles & BWTLocateBitmap)
  {
    unsigned long nextLocate = pos;
    unsigned locateOffset = 0;
    while (!gt_BWTSeqPosHasLocateInfo(bwtSeq, nextLocate, extBits))
      nextLocate = BWTSeqLFMap(bwtSeq, nextLocate, extBits), ++locateOffset;
    EISRetrieveExtraBits(bwtSeq->seqIdx, nextLocate,
                         EBRF_RETRIEVE_CWBITS | EBRF_RETRIEVE_VARBITS,
                         extBits, bwtSeq->hint);
    {
      unsigned long maxPosVal = ((bwtSeq->featureToggles & BWTReversiblySorted)?
                          (BWTSeqLength(bwtSeq) - 1)
                          /bwtSeq->locateSampleInterval:
                          BWTSeqLength(bwtSeq) - 1);
      unsigned bitsPerCount = requiredUlongBits(extBits->len),
        bitsPerBWTPos = requiredUlongBits(extBits->len - 1),
        bitsPerOrigPos = requiredUlongBits(maxPosVal);
      BitOffset locateRecordIndex =
        gt_bs1BitsCount(extBits->cwPart, extBits->cwOffset,
                     nextLocate - extBits->start),
        locateRecordOffset = ((bwtSeq->featureToggles & BWTLocateCount?
                               bitsPerBWTPos:0) + bitsPerOrigPos)
        * locateRecordIndex
        + ((bwtSeq->featureToggles & BWTLocateCount)?bitsPerCount:0);
      unsigned long matchPos =
        gt_bsGetUlong(
          extBits->varPart, extBits->varOffset + locateRecordOffset
          + ((bwtSeq->featureToggles & BWTLocateCount)?bitsPerBWTPos:0),
          bitsPerOrigPos);
      if (bwtSeq->featureToggles & BWTReversiblySorted)
        matchPos = matchPos * bwtSeq->locateSampleInterval;
      matchPos += locateOffset;
      gt_assert(!(bwtSeq->featureToggles & BWTLocateCount)
             || gt_bsGetUlong(extBits->varPart,
                            extBits->varOffset + locateRecordOffset,
                            bitsPerBWTPos)
             == nextLocate - extBits->start);
      return matchPos;
    }
  }
  else if (bwtSeq->featureToggles & BWTLocateCount)
  {
    BitOffset markOffset;
    unsigned long nextLocate = pos, matchPos;
    unsigned locateOffset = 0;  /* mark is at most locateInterval
                                 * positions away */
    unsigned bitsPerOrigPos
      = requiredUlongBits(((bwtSeq->featureToggles & BWTReversiblySorted)?
                            (BWTSeqLength(bwtSeq) - 1)
                            /bwtSeq->locateSampleInterval:
                            BWTSeqLength(bwtSeq) - 1));
    while ((markOffset = searchLocateCountMark(bwtSeq, nextLocate,
                                               extBits)) == 0)
    {
      nextLocate = BWTSeqLFMap(bwtSeq, nextLocate, extBits);
      ++locateOffset;
      gt_assert(locateOffset <= BWTSeqLength(bwtSeq));
    }
    matchPos = gt_bsGetUlong(extBits->varPart, markOffset, bitsPerOrigPos);
    if (bwtSeq->featureToggles & BWTReversiblySorted)
        matchPos = matchPos * bwtSeq->locateSampleInterval;
    matchPos += locateOffset;
    return matchPos;
  }
  /* Internal error: Trying to locate in BWT sequence index without locate
     information. */
   abort();
   return 0; /* shut up compiler */
}

static inline BitOffset
locateVarBits(const BWTSeq *bwtSeq, struct extBitsRetrieval *extBits)
{
  BitOffset numLocBits = 0;
  unsigned bitsPerBWTPos = requiredUlongBits(extBits->len - 1),
    bitsPerOrigPos = requiredUlongBits(
      ((bwtSeq->featureToggles & BWTReversiblySorted) ?
       (BWTSeqLength(bwtSeq) - 1) / bwtSeq->locateSampleInterval :
       BWTSeqLength(bwtSeq) - 1));
  if (bwtSeq->featureToggles & BWTLocateBitmap)
  {
    unsigned numMarks = gt_bs1BitsCount(extBits->cwPart, extBits->cwOffset,
                                     extBits->len);
    numLocBits = numMarks * bitsPerOrigPos;
  }
  else if (bwtSeq->featureToggles & BWTLocateCount)
  {
    BitOffset markOffset = extBits->varOffset;
    unsigned bitsPerCount = requiredUlongBits(extBits->len);
    unsigned numMarks = gt_bsGetUlong(extBits->varPart, markOffset,
                                    bitsPerCount);
    numLocBits = bitsPerCount + numMarks * (bitsPerBWTPos + bitsPerOrigPos);
  }
  return numLocBits;
}

unsigned long
gt_BWTSeqGetRankSort(const BWTSeq *bwtSeq, unsigned long pos,
                  AlphabetRangeID range, struct extBitsRetrieval *extBits)
{
  BitOffset locVarBits;
  AlphabetRangeSize rSize;
  gt_assert(bwtSeq->rangeSort[range] == SORTMODE_RANK);
  EISRetrieveExtraBits(bwtSeq->seqIdx, pos, EBRF_RETRIEVE_CWBITS
                       | EBRF_RETRIEVE_VARBITS, extBits, bwtSeq->hint);
  locVarBits = locateVarBits(bwtSeq, extBits);
  rSize = MRAEncGetRangeSize(EISGetAlphabet(bwtSeq->seqIdx), range);
  {
    unsigned long ranks[rSize * 2];
    unsigned long BWTRankTotal = 0;
    AlphabetRangeSize i;
    EISPosPairRangeRank(bwtSeq->seqIdx, range, extBits->start, pos, ranks,
                        bwtSeq->hint);
    for (i = 0; i < rSize; ++i)
      BWTRankTotal += ranks[i + rSize] - ranks[i];
    return gt_bsGetUlong(extBits->varPart, extBits->varOffset + locVarBits
                       + bwtSeq->bitsPerOrigRank * BWTRankTotal,
                       bwtSeq->bitsPerOrigRank);
  }
}

void
gt_BWTSeqInitLocateHandling(BWTSeq *bwtSeq,
                         const enum rangeSortMode *defaultRangeSort)
{
  struct encIdxSeq *seqIdx;
  struct locateHeader locHeader;
  gt_assert(bwtSeq);
  seqIdx = bwtSeq->seqIdx;
  if (!readLocateInfoHeader(seqIdx, &locHeader)
      || !locHeader.locateInterval)
  {
    gt_log_log("Index does not contain locate information.\n"
            "Localization of matches will not be supported!");
    bwtSeq->locateSampleInterval = 0;
    bwtSeq->featureToggles = BWTBaseFeatures;
  }
  else
  {
    bwtSeq->locateSampleInterval = locHeader.locateInterval;
    bwtSeq->rot0Pos = locHeader.rot0Pos;
    /* FIXME: this really deserves its own header */
    bwtSeq->featureToggles = locHeader.featureToggles;

    if (readRankSortHeader(seqIdx, &bwtSeq->bitsPerOrigRank,
                           bwtSeq->alphabet, bwtSeq->rangeSort))
      ;
    else
    {
      AlphabetRangeID numRanges = MRAEncGetNumRanges(bwtSeq->alphabet);
      bwtSeq->bitsPerOrigRank = 0;
      memcpy(bwtSeq->rangeSort, defaultRangeSort,
             numRanges * sizeof (defaultRangeSort[0]));
    }
  }
}
