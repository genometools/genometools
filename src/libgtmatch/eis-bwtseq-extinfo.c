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

#include "libgtcore/chardef.h"
#include "libgtcore/log.h"
#include "libgtcore/minmax.h"
#include "libgtcore/unused.h"
#include "libgtmatch/encseq-specialsrank.h"
#include "libgtmatch/eis-bitpackseqpos.h"
#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseq-extinfo.h"
#include "libgtmatch/eis-bwtseq-priv.h"
#include "libgtmatch/eis-bwtseq-context.h"
#include "libgtmatch/eis-headerid.h"
#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-sa-common.h"

/**
 * @file eis-bwtseq-extinfo.c generic methods for bwt index creation and
 * everything related to storing/retrieving locate information
 */

struct locateHeader
{
  Seqpos rot0Pos;
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
  DefinedSeqpos rot0Pos;
  const struct locateHeaderWriteInfo *headerSrc = cbData;
  assert(cbData);
  headerData.locateInterval = headerSrc->locateInterval;
  rot0Pos = SASSGetRot0Pos(headerSrc->src);
  if (!rot0Pos.defined)
  {
    fputs("Invalid index construction: position of suffix 0 unknown!\n",
          stderr);
    abort();
  }
  headerData.rot0Pos = rot0Pos.valueseqpos;
  headerData.featureToggles = headerSrc->featureToggles;
  return fwrite(&headerData, sizeof (headerData), 1, fp);
}

/**
 * @return 0 on error
 */
static inline int
readLocateInfoHeader(EISeq *seqIdx, struct locateHeader *headerData)
{
  FILE *fp;
  assert(seqIdx && headerData);
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
  assert(cbData);
  if (fwrite(&headerData->bitsPerOrigRank, sizeof (headerData->bitsPerOrigRank),
             1, fp) != 1)
    return 0;
  {
    size_t i, numRanges = MRAEncGetNumRanges(headerData->alphabet);
    for (i = 0; i < numRanges; ++i)
    {
      int16_t mode = headerData->rangeSort[i];
      if (fwrite(&mode, sizeof (mode), 1, fp) != 1)
        return 0;
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
  assert(seqIdx && alphabet && bitsPerOrigRank && rangeSort);
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
  Seqpos bwtPos, origPos;
};

struct addLocateInfoState
{
  Seqpos seqLen, extraLocMarksUpperBound;
  const MRAEnc *alphabet;
  const enum rangeSortMode *rangeSort;
  SeqDataReader readSeqpos;
  RandomSeqAccessor origSeqAccess;
  unsigned locateInterval, bitsPerOrigPos, bitsPerSeqpos,
    bitsPerOrigRank;
  const SpecialsRankLookup *sprTable;
  int featureToggles;
  size_t revMapQueueSize;
  struct seqRevMapEntry *revMapQueue;
  size_t origRanksQueueSize;
  Seqpos *origRanksQueue;
  BWTSeqContextRetrieverFactory *ctxFactory;
};

static inline unsigned
bitsPerPosUpperBoundWithoutSegmentAggregates(struct addLocateInfoState *state)
{
  unsigned bitsPerPos = 0;
  if (state->featureToggles & BWTLocateCount)
    /* bitsPerSeqpos is actually an upper bound of log(segmentLen) to
     * identify each locate mark */
    bitsPerPos += state->bitsPerSeqpos;
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
  assert(cbState);
  if (state->featureToggles & (BWTLocateCount | BWTLocateBitmap))
  {
    unsigned maxBitsPerPos = bitsPerPosUpperBound(state);
    result->maxBitsPerPos = maxBitsPerPos;
    if (desc)
    {
      size_t i, maxSegLen = 0;
      BitOffset maxBitsTotal = 0;
      Seqpos numSegmentsTotal = 0;
      unsigned bitsPerOrigRank = state->bitsPerOrigRank;
      for (i = 0; i < numSegmentDesc; ++i)
      {
        size_t len = desc[i].len;
        maxSegLen = MAX(len, maxSegLen);
        if (state->featureToggles & BWTLocateCount)
          maxBitsTotal += requiredSeqposBits(len) * desc[i].repeatCount;
        numSegmentsTotal += desc[i].repeatCount;
      }
      maxBitsTotal +=
        (state->seqLen / state->locateInterval + state->extraLocMarksUpperBound)
        * (((state->featureToggles & BWTLocateCount) ?
            requiredSeqposBits(maxSegLen) : 0)
           + state->bitsPerOrigPos);
      if (bitsPerOrigRank)
      {
        Seqpos maxRank = specialsRank(state->sprTable, state->seqLen);
        /* rank values stored */
        maxBitsTotal += maxRank * bitsPerOrigRank;
      }
      result->maxBitsPerBucket =
        /* locate marks and rank sorts */
        maxSegLen * bitsPerPosUpperBoundWithoutSegmentAggregates(state)
        /* count of locate marks */
        + ((state->featureToggles & BWTLocateCount) ?
           requiredSeqposBits(maxSegLen) : 0);
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
                       SeqDataReader readSeqpos,
                       const MRAEnc *alphabet, const struct seqStats *stats,
                       const enum rangeSortMode *rangeSort,
                       Seqpos srcLen, const struct bwtParam *params,
                       const SpecialsRankLookup *sprTable,
                       unsigned bitsPerOrigRank,
                       BWTSeqContextRetrieverFactory *ctxFactory)
{
  Seqpos lastPos;
  unsigned aggregationExpVal;
  unsigned locateInterval;
  assert(state);
  state->alphabet = alphabet;
  state->seqLen = srcLen;
  state->rangeSort = rangeSort;
  state->readSeqpos = readSeqpos;
  state->origSeqAccess = origSeqAccess;
  state->featureToggles = params->featureToggles;
  aggregationExpVal = estimateSegmentSize(&params->seqParams);
  locateInterval = params->locateInterval;
  lastPos = srcLen - 1;
  state->locateInterval = locateInterval;
  state->bitsPerSeqpos = requiredSeqposBits(lastPos);
  state->bitsPerOrigRank = bitsPerOrigRank;
  state->sprTable = sprTable;
  if (locateInterval)
  {
    if (params->featureToggles & BWTReversiblySorted)
      state->bitsPerOrigPos = requiredSeqposBits(lastPos/locateInterval);
    else
      state->bitsPerOrigPos = requiredSeqposBits(lastPos);
    state->origRanksQueueSize = state->revMapQueueSize = aggregationExpVal;
    state->revMapQueue = ma_malloc(sizeof (state->revMapQueue[0])
                                   * state->revMapQueueSize);
    state->origRanksQueue = ma_malloc(sizeof (state->origRanksQueue[0])
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
      Seqpos stdLocMarks = srcLen / locateInterval,
        extraLocMarksUpperBound = MIN(srcLen/2, srcLen - stdLocMarks);
      if (stats)
      {
        Seqpos nonValSortSyms = 0;
        unsigned i;
        for (i = 0; i <= UINT8_MAX; ++i)
          if (MRAEncSymbolHasValidMapping(alphabet, i)
              && !MRAEncSymbolIsInSelectedRanges(
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
  ma_free(state->revMapQueue);
  ma_free(state->origRanksQueue);
}

static inline int
isSortModeTransition(RandomSeqAccessor origSeqAccess, Seqpos seqLen,
                     const MRAEnc *alphabet,
                     const enum rangeSortMode *rangeSort,
                     Seqpos pos)
{
  Symbol syms[2];
  if (pos > 0 && pos < seqLen - 1)
  {
#ifndef NDEBUG
    int retcode =
#endif
      accessSequence(origSeqAccess, syms, pos - 1, 2);
    assert(retcode == 2);
  }
  else if (pos == 0)
  {
#ifndef NDEBUG
    int retcode =
#endif
      accessSequence(origSeqAccess, syms + 1, 0, 1);
    assert(retcode == 1);
    syms[0] = UNDEFBWTCHAR;
  }
  else /* pos == seqLen - 1 */
  {
#ifndef NDEBUG
    int retcode =
#endif
      accessSequence(origSeqAccess, syms, seqLen - 2, 1);
    assert(retcode == 1);
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
              UNUSED Seqpos start, Seqpos len, void *cbState, Error *err)
{
  BitOffset bitsWritten = 0;
  struct addLocateInfoState *state = cbState;
  unsigned bitsPerBWTPos, bitsPerOrigPos, bitsPerOrigRank;
  int retcode;
  unsigned locateInterval;
  assert(varDest && cbState);
  locateInterval = state->locateInterval;
  /* 0. resize caches if necessary */
  if (locateInterval && len > state->revMapQueueSize)
  {
    state->revMapQueue = ma_realloc(state->revMapQueue,
                                    len * sizeof (state->revMapQueue[0]));
    state->revMapQueueSize = len;
  }
  if (state->sprTable && len > state->origRanksQueueSize)
  {
    state->origRanksQueue = ma_realloc(state->origRanksQueue,
                                       len * sizeof (state->origRanksQueue[0]));
    state->origRanksQueueSize = len;
  }
  bitsPerBWTPos = requiredSeqposBits(len - 1);
  bitsPerOrigPos = state->bitsPerOrigPos;
  bitsPerOrigRank = state->bitsPerOrigRank;
  {
    Seqpos i, mapVal;
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
        if ((retcode = SDRRead(state->readSeqpos, &mapVal, 1, err)) != 1)
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
            bsClearBit(cwDest, cwOffset + i);
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
        BWTSCRFMapAdvance(state->ctxFactory, &mapVal, 1);
      }
    }
    /* 2. copy revMapQueue into output */
    if (locateInterval)
    {
      int locateCount = state->featureToggles & BWTLocateCount;
      if (locateCount)
      {
        unsigned bitsPerCount = requiredSeqposBits(len);
        bsStoreSeqpos(varDest, varOffset + bitsWritten, bitsPerCount,
                      revMapQueueLen);
        bitsWritten += bitsPerCount;
      }
      for (i = 0; i < revMapQueueLen; ++i)
      {
        if (locateCount)
        {
          bsStoreSeqpos(varDest, varOffset + bitsWritten, bitsPerBWTPos,
                        state->revMapQueue[i].bwtPos);
          bitsWritten += bitsPerBWTPos;
        }
        bsStoreSeqpos(varDest, varOffset + bitsWritten, bitsPerOrigPos,
                      state->revMapQueue[i].origPos);
        bitsWritten += bitsPerOrigPos;
      }
    }
    if (bitsPerOrigRank)
    {
      unsigned bitsPerOrigRank = state->bitsPerOrigRank;
      if (origRanksQueueLen)
      {
        bsStoreUniformSeqposArray(varDest, varOffset + bitsWritten,
                                  bitsPerOrigRank, origRanksQueueLen,
                                  state->origRanksQueue);
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

extern EISeq *
createBWTSeqGeneric(const struct bwtParam *params, indexCreateFunc createIndex,
                    SASeqSrc *src,
                    const enum rangeSortMode rangeSort[],
                    const SpecialsRankLookup *sprTable,
                    Error *err)
{
  struct encIdxSeq *baseSeqIdx = NULL;
  struct addLocateInfoState varState;
  bool varStateIsInitialized = false;
  unsigned locateInterval;
  BWTSeqContextRetrieverFactory *buildContextMap = NULL;
  assert(src && params && err);
  error_check(err);
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
    Seqpos totalLen = SASSGetLength(src);
    const MRAEnc *alphabet = SASSGetMRAEnc(src);
    MRAEnc *baseAlphabet = SASSNewMRAEnc(src);
    /* FIXME: this  has to work also when locateInterval == 0 and
     * sprTable != NULL */
    if (params->ctxMapILog != CTX_MAP_ILOG_NOMAP)
      buildContextMap = newBWTSeqContextRetrieverFactory(totalLen,
                                                         params->ctxMapILog);
    if (locateInterval)
    {
      ++numHeaders;
      if (sortModeHeaderNeeded(alphabet, rangeSort, sprTable))
      {
        Seqpos
#ifndef NDEBUG
          origSeqLen = getencseqtotallength(SPRTGetOrigEncseq(sprTable)),
#endif
          maxRank;
        assert(origSeqLen == totalLen - 1);
        maxRank = specialsRank(sprTable, totalLen - 1);
        bitsPerOrigRank = sortModeHeader.bitsPerOrigRank
          = requiredSeqposBits(maxRank);
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
          error_set(err, "error: locate sampling requested but not available"
                    " for project %s\n", str_get(params->projectName));
        }
      }
    }
    if (!(baseSeqIdx
          = createIndex(totalLen, params->projectName, baseAlphabet,
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
      if (!BWTSCRFFinished(buildContextMap))
      {
        fputs("error: context table construction incomplete!\n", stderr);
      }
      else
      {
        BWTSeqContextRetriever *ctxRetrieve =
          BWTSCRFGet(buildContextMap, NULL, params->projectName);
        deleteBWTSeqCR(ctxRetrieve);
      }
    }
  } while (0);
  if (buildContextMap)
    deleteBWTSeqContextRetrieverFactory(buildContextMap);
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
      bitsPerOrigPos = requiredSeqposBits(
        ((bwtSeq->featureToggles & BWTReversiblySorted) ?
         (BWTSeqLength(bwtSeq) - 1) / bwtSeq->locateSampleInterval :
         BWTSeqLength(bwtSeq) - 1));
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
      nextLocate = BWTSeqLFMap(bwtSeq, nextLocate, extBits), ++locateOffset;
    EISRetrieveExtraBits(bwtSeq->seqIdx, nextLocate,
                         EBRF_RETRIEVE_CWBITS | EBRF_RETRIEVE_VARBITS,
                         extBits, bwtSeq->hint);
    {
      Seqpos maxPosVal = ((bwtSeq->featureToggles & BWTReversiblySorted)?
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
      if (bwtSeq->featureToggles & BWTReversiblySorted)
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
      = requiredSeqposBits(((bwtSeq->featureToggles & BWTReversiblySorted)?
                            (BWTSeqLength(bwtSeq) - 1)
                            /bwtSeq->locateSampleInterval:
                            BWTSeqLength(bwtSeq) - 1));
    while ((markOffset = searchLocateCountMark(bwtSeq, nextLocate,
                                               extBits)) == 0)
    {
      nextLocate = BWTSeqLFMap(bwtSeq, nextLocate, extBits);
      ++locateOffset;
    }
    matchPos = bsGetSeqpos(extBits->varPart, markOffset, bitsPerOrigPos);
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
  unsigned bitsPerBWTPos = requiredSeqposBits(extBits->len - 1),
    bitsPerOrigPos = requiredSeqposBits(
      ((bwtSeq->featureToggles & BWTReversiblySorted) ?
       (BWTSeqLength(bwtSeq) - 1) / bwtSeq->locateSampleInterval :
       BWTSeqLength(bwtSeq) - 1));
  if (bwtSeq->featureToggles & BWTLocateBitmap)
  {
    unsigned numMarks = bs1BitsCount(extBits->cwPart, extBits->cwOffset,
                                     extBits->len);
    numLocBits = numMarks * bitsPerOrigPos;
  }
  else if (bwtSeq->featureToggles & BWTLocateCount)
  {
    BitOffset markOffset = extBits->varOffset;
    unsigned bitsPerCount = requiredSeqposBits(extBits->len);
    unsigned numMarks = bsGetSeqpos(extBits->varPart, markOffset,
                                    bitsPerCount);
    numLocBits = bitsPerCount + numMarks * (bitsPerBWTPos + bitsPerOrigPos);
  }
  return numLocBits;
}

extern Seqpos
BWTSeqGetRankSort(const BWTSeq *bwtSeq, Seqpos pos, AlphabetRangeID range,
                  struct extBitsRetrieval *extBits)
{
  BitOffset locVarBits;
  AlphabetRangeSize rSize;
  assert(bwtSeq->rangeSort[range] == SORTMODE_RANK);
  EISRetrieveExtraBits(bwtSeq->seqIdx, pos, EBRF_RETRIEVE_CWBITS
                       | EBRF_RETRIEVE_VARBITS, extBits, bwtSeq->hint);
  locVarBits = locateVarBits(bwtSeq, extBits);
  rSize = MRAEncGetRangeSize(EISGetAlphabet(bwtSeq->seqIdx), range);
  {
    Seqpos ranks[rSize * 2];
    Seqpos BWTRankTotal = 0;
    AlphabetRangeSize i;
    EISPosPairRangeRank(bwtSeq->seqIdx, range, extBits->start, pos, ranks,
                        bwtSeq->hint);
    for (i = 0; i < rSize; ++i)
      BWTRankTotal += ranks[i + rSize] - ranks[i];
    return bsGetSeqpos(extBits->varPart, extBits->varOffset + locVarBits
                       + bwtSeq->bitsPerOrigRank * BWTRankTotal,
                       bwtSeq->bitsPerOrigRank);
  }
}

extern void
BWTSeqInitLocateHandling(BWTSeq *bwtSeq,
                         const enum rangeSortMode *defaultRangeSort)
{
  struct encIdxSeq *seqIdx;
  struct locateHeader locHeader;
  assert(bwtSeq);
  seqIdx = bwtSeq->seqIdx;
  if (!readLocateInfoHeader(seqIdx, &locHeader)
      || !locHeader.locateInterval)
  {
    log_log("Index does not contain locate information.\n"
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
