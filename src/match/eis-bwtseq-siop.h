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

#ifndef EIS_BWTSEQ_SIOP_H
#define EIS_BWTSEQ_SIOP_H

/* trivial operations on BWTSeq objects go here for speed */

#include "match/eis-bwtseq.h"
#include "match/eis-bwtseq-priv.h"
#include "stamp.h"

static inline const MRAEnc *
BWTSeqGetAlphabet(const BWTSeq *seq)
{
  return seq->alphabet;
}

static inline unsigned long
BWTSeqLength(const BWTSeq *seq)
{
  return EISLength(seq->seqIdx);
}

static inline unsigned long
BWTSeqTerminatorPos(const BWTSeq *bwtSeq)
{
  return bwtSeq->rot0Pos;
}

static inline bool
BWTSeqHasLocateInformation(const BWTSeq *bwtSeq)
{
  return bwtSeq->locateSampleInterval != 0;
}

static inline unsigned long
BWTSeqTransformedOcc(const BWTSeq *bwtSeq, Symbol tsym, unsigned long pos)
{
  gt_assert(bwtSeq);
  /* two counts must be treated specially:
   * 1. for the symbols mapped to the same value as the terminator
   * 2. for queries of the terminator itself */
  if (tsym < bwtSeq->bwtTerminatorFallback)
    return EISSymTransformedRank(bwtSeq->seqIdx, tsym, pos, bwtSeq->hint);
  else if (tsym > bwtSeq->bwtTerminatorFallback
           && tsym != bwtSeq->alphabetSize - 1)
    return EISSymTransformedRank(bwtSeq->seqIdx, tsym, pos, bwtSeq->hint);
  else if (tsym == bwtSeq->bwtTerminatorFallback)
    return EISSymTransformedRank(bwtSeq->seqIdx, tsym, pos, bwtSeq->hint)
      - ((pos > BWTSeqTerminatorPos(bwtSeq))?1:0);
  else /* tsym == not flattened terminator == alphabetSize - 1 */
  {
    gt_assert(tsym == bwtSeq->alphabetSize - 1);
    return (pos > BWTSeqTerminatorPos(bwtSeq))?1:0;
  }
}

static inline GtUlongPair
BWTSeqTransformedPosPairOcc(const BWTSeq *bwtSeq, Symbol tSym,
                            unsigned long posA, unsigned long posB)
{
  gt_assert(bwtSeq);
  /* two counts must be treated specially:
   * 1. for the symbols mapped to the same value as the terminator
   * 2. for queries of the terminator itself */
  if (tSym < bwtSeq->bwtTerminatorFallback)
    return EISSymTransformedPosPairRank(bwtSeq->seqIdx, tSym, posA, posB,
                                        bwtSeq->hint);
  else if (tSym > bwtSeq->bwtTerminatorFallback
           && tSym != bwtSeq->alphabetSize - 1)
    return EISSymTransformedPosPairRank(bwtSeq->seqIdx, tSym, posA, posB,
                                        bwtSeq->hint);
  else if (tSym == bwtSeq->bwtTerminatorFallback)
    return EISSymTransformedPosPairRank(bwtSeq->seqIdx, tSym, posA, posB,
                                        bwtSeq->hint);
/*       - ((pos > BWTSeqTerminatorPos(bwtSeq))?1:0); */
  else /* tSym == not flattened terminator == alphabetSize - 1 */
  {
    GtUlongPair occ;
    gt_assert(tSym == bwtSeq->alphabetSize - 1);
    occ.a = (posA > BWTSeqTerminatorPos(bwtSeq))?1:0;
    occ.b = (posB > BWTSeqTerminatorPos(bwtSeq))?1:0;
    return occ;
  }
}

static inline unsigned long
BWTSeqOcc(const BWTSeq *bwtSeq, Symbol sym, unsigned long pos)
{
  gt_assert(bwtSeq);
  Symbol tSym = MRAEncMapSymbol(BWTSeqGetAlphabet(bwtSeq), sym);
  return BWTSeqTransformedOcc(bwtSeq, tSym, pos);
}

static inline GtUlongPair
BWTSeqPosPairOcc(const BWTSeq *bwtSeq, Symbol sym, unsigned long posA,
                 unsigned long posB)
{
  gt_assert(bwtSeq);
  Symbol tSym = MRAEncMapSymbol(BWTSeqGetAlphabet(bwtSeq), sym);
  return BWTSeqTransformedPosPairOcc(bwtSeq, tSym, posA, posB);
}

static inline void
BWTSeqRangeOcc(const BWTSeq *bwtSeq, AlphabetRangeID range, unsigned long pos,
               unsigned long *rangeOccs)
{
  gt_assert(bwtSeq && rangeOccs);
  gt_assert(range < MRAEncGetNumRanges(BWTSeqGetAlphabet(bwtSeq)));
  EISRangeRank(bwtSeq->seqIdx, range, pos, rangeOccs, bwtSeq->hint);
  if (range == bwtSeq->bwtTerminatorFallbackRange)
  {
    const MRAEnc *alphabet = BWTSeqGetAlphabet(bwtSeq);
    AlphabetRangeSize
      termIdx = MRAEncMapSymbol(alphabet, bwtTerminatorSym)
      - MRAEncGetRangeBase(alphabet, range),
      fbIdx = bwtSeq->bwtTerminatorFallback
      - MRAEncGetRangeBase(alphabet, range);
    if (pos > BWTSeqTerminatorPos(bwtSeq))
    {
      rangeOccs[termIdx] = 1;
      --(rangeOccs[fbIdx]);
    }
    else
      rangeOccs[termIdx] = 0;
  }
}

static inline void
BWTSeqPosPairRangeOcc(const BWTSeq *bwtSeq, AlphabetRangeID range,
                      unsigned long posA, unsigned long posB,
                      unsigned long *rangeOccs)
{
  gt_assert(bwtSeq && rangeOccs);
  gt_assert(posA <= posB);
  gt_assert(range < MRAEncGetNumRanges(BWTSeqGetAlphabet(bwtSeq)));
  EISPosPairRangeRank(bwtSeq->seqIdx, range, posA, posB, rangeOccs,
                      bwtSeq->hint);
  if (range == bwtSeq->bwtTerminatorFallbackRange)
  {
    const MRAEnc *alphabet = BWTSeqGetAlphabet(bwtSeq);
    AlphabetRangeSize rSize = MRAEncGetRangeSize(alphabet, range),
      termIdx = MRAEncMapSymbol(alphabet, bwtTerminatorSym)
      - MRAEncGetRangeBase(alphabet, range),
      fbIdx = bwtSeq->bwtTerminatorFallback
      - MRAEncGetRangeBase(alphabet, range);
    memmove(rangeOccs + termIdx + 1, rangeOccs + termIdx,
            sizeof (rangeOccs[0]) * (rSize - 1));
    rangeOccs[termIdx] = rangeOccs[termIdx + rSize] = 0;
    if (posB > BWTSeqTerminatorPos(bwtSeq))
    {
      rangeOccs[termIdx + rSize] = 1;
      --(rangeOccs[fbIdx + rSize]);
      if (posA > BWTSeqTerminatorPos(bwtSeq))
      {
        rangeOccs[termIdx] = 1;
        --(rangeOccs[fbIdx]);
      }
      else
        rangeOccs[termIdx] = 0;
    }
    else
      rangeOccs[rSize + termIdx] = rangeOccs[termIdx] = 0;
  }
}

static inline unsigned long
BWTSeqLFMap(const BWTSeq *bwtSeq, unsigned long LPos,
            struct extBitsRetrieval *extBits)
{
  Symbol tSym = EISGetTransformedSym(bwtSeq->seqIdx, LPos, bwtSeq->hint);
  unsigned long FPos;
  const MRAEnc *alphabet = BWTSeqGetAlphabet(bwtSeq);
  if (LPos != BWTSeqTerminatorPos(bwtSeq))
  {
    AlphabetRangeID range = MRAEncGetRangeOfSymbol(alphabet, tSym);
    switch (bwtSeq->rangeSort[range])
    {
    case SORTMODE_VALUE:
      FPos = bwtSeq->count[tSym] + BWTSeqTransformedOcc(bwtSeq, tSym, LPos);
      break;
    case SORTMODE_RANK:
      FPos = bwtSeq->count[MRAEncGetRangeBase(alphabet, range)]
        + gt_BWTSeqGetRankSort(bwtSeq, LPos, range, extBits);
      break;
    case SORTMODE_UNDEFINED:
    default:
#ifndef _NDEBUG
      fputs("Requesting LF-map for symbol from range of undefined sorting.\n",
            stderr);
      abort();
#endif
      FPos = (unsigned long)-1;
      break;
    }
  }
  else
  {
    gt_assert(tSym == bwtSeq->bwtTerminatorFallback);
    FPos = bwtSeq->count[bwtSeq->alphabetSize - 1];
  }
  return FPos;
}

static inline unsigned long
BWTSeqAggCount(const BWTSeq *bwtSeq, Symbol sym)
{
  Symbol  tSym;
  gt_assert(bwtSeq);
  gt_assert(MRAEncSymbolHasValidMapping(BWTSeqGetAlphabet(bwtSeq), sym));
  tSym = MRAEncMapSymbol(BWTSeqGetAlphabet(bwtSeq), sym);
  return bwtSeq->count[tSym];
}

static inline unsigned long
BWTSeqAggTransformedCount(const BWTSeq *bwtSeq, Symbol tSym)
{
  gt_assert(bwtSeq);
  gt_assert(tSym <= bwtSeq->alphabetSize);
  return bwtSeq->count[tSym];
}

static inline struct matchBound *
BWTSeqIncrMatch(const BWTSeq *bwtSeq, struct matchBound *limits,
                Symbol nextSym)
{
  const MRAEnc *alphabet;
  Symbol curSym;
  gt_assert(bwtSeq && limits);
  gt_assert(limits->start < bwtSeq->count[bwtSeq->alphabetSize]
         && limits->end < bwtSeq->count[bwtSeq->alphabetSize]);
  alphabet = BWTSeqGetAlphabet(bwtSeq);
  curSym = MRAEncMapSymbol(alphabet, nextSym);
  gt_assert(MRAEncSymbolHasValidMapping(alphabet, curSym));
  limits->start = bwtSeq->count[curSym]
    + BWTSeqTransformedOcc(bwtSeq, curSym, limits->start);
  limits->end = bwtSeq->count[curSym]
    + BWTSeqTransformedOcc(bwtSeq, curSym, limits->end);
  return limits;
}

static inline const EISeq *
BWTSeqGetEncIdxSeq(const BWTSeq *bwtSeq)
{
  gt_assert(bwtSeq);
  return bwtSeq->seqIdx;
}

static inline Symbol
BWTSeqGetSym(const BWTSeq *bwtSeq, unsigned long pos)
{
  gt_assert(bwtSeq);
  return EISGetSym(bwtSeq->seqIdx, pos, bwtSeq->hint);
}

static inline bool
EMIGetNextMatch(struct BWTSeqExactMatchesIterator *iter, unsigned long *pos,
                const BWTSeq *bwtSeq)
{
  if (iter->nextMatchBWTPos < iter->bounds.end)
  {
    /*printf("nextMatchBWTPos=%lu\n",(unsigned long) iter->nextMatchBWTPos);*/
    *pos = gt_BWTSeqLocateMatch(bwtSeq, iter->nextMatchBWTPos, &iter->extBits);
    iter->nextMatchBWTPos++;
    return true;
  }
  else
    return false;
}

#endif
