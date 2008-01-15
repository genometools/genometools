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
#ifndef EIS_BWTSEQSIMPLEOP_H
#define EIS_BWTSEQSIMPLEOP_H

#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseqpriv.h"

static inline const MRAEnc *
BWTSeqGetAlphabet(const BWTSeq *seq)
{
  return seq->alphabet;
}

static inline Seqpos
BWTSeqLength(const BWTSeq *seq)
{
  return EISLength(seq->seqIdx);
}

static inline bool
BWTSeqHasLocateInformation(const BWTSeq *bwtSeq)
{
  return bwtSeq->locateSampleInterval != 0;
}

static inline Seqpos
BWTSeqTransformedOcc(const BWTSeq *bwtSeq, Symbol tsym, Seqpos pos)
{
  assert(bwtSeq);
  /* two counts must be treated specially:
   * 1. for the symbols mapped to the same value as the terminator
   * 2. for queries of the terminator itself */
  if (tsym < bwtSeq->bwtTerminatorFallback)
    return EISSymTransformedRank(bwtSeq->seqIdx, tsym, pos, bwtSeq->hint);
  else if (tsym > bwtSeq->bwtTerminatorFallback
           && tsym != bwtSeq->alphabetSize - 1)
    return EISSymTransformedRank(bwtSeq->seqIdx, tsym, pos, bwtSeq->hint);
  else if (tsym == bwtSeq->bwtTerminatorFallback)
    return EISSymTransformedRank(bwtSeq->seqIdx, tsym, pos, bwtSeq->hint);
/*       - ((pos > bwtSeq->longest)?1:0); */
  else /* tsym == not flattened terminator == alphabetSize - 1 */
  {
    assert(tsym == bwtSeq->alphabetSize - 1);
    return (pos > bwtSeq->longest)?1:0;
  }
}

static inline Seqpos
BWTSeqOcc(const BWTSeq *bwtSeq, Symbol sym, Seqpos pos)
{
  assert(bwtSeq);
  Symbol tsym = MRAEncMapSymbol(bwtSeq->alphabet, sym);
  return BWTSeqTransformedOcc(bwtSeq, tsym, pos);
}

static inline Seqpos
BWTSeqLFMap(const BWTSeq *bwtSeq, Seqpos pos)
{
  Symbol tSym = EISGetTransformedSym(bwtSeq->seqIdx, pos, bwtSeq->hint);
  if (pos != bwtSeq->longest)
  {
    return bwtSeq->count[tSym] + BWTSeqTransformedOcc(bwtSeq, tSym, pos);
  }
  else
  {
    assert(tSym == bwtSeq->bwtTerminatorFallback);
    return bwtSeq->count[bwtSeq->alphabetSize - 1];
  }
}

static inline Seqpos
BWTSeqAggCount(const BWTSeq *bwtSeq, Symbol sym)
{
  Symbol  tSym;
  assert(bwtSeq);
  assert(MRAEncSymbolHasValidMapping(BWTSeqGetAlphabet(bwtSeq), sym));
  tSym = MRAEncMapSymbol(BWTSeqGetAlphabet(bwtSeq), sym);
  return bwtSeq->count[tSym];
}

static inline Seqpos
BWTSeqAggTransformedCount(const BWTSeq *bwtSeq, Symbol tSym)
{
  assert(bwtSeq);
  assert(tSym <= bwtSeq->alphabetSize);
  return bwtSeq->count[tSym];
}

static struct matchBound *
BWTSeqIncrMatch(const BWTSeq *bwtSeq, struct matchBound *limits,
                Symbol nextSym)
{
  const MRAEnc *alphabet;
  Symbol curSym;
  assert(bwtSeq && limits);
  assert(limits->upper < bwtSeq->count[bwtSeq->alphabetSize]
         && limits->lower < bwtSeq->count[bwtSeq->alphabetSize]);
  alphabet = BWTSeqGetAlphabet(bwtSeq);
  curSym = MRAEncMapSymbol(alphabet, nextSym);
  assert(MRAEncSymbolHasValidMapping(alphabet, curSym));
  limits->upper = bwtSeq->count[curSym]
    + BWTSeqTransformedOcc(bwtSeq, curSym, limits->upper);
  limits->lower = bwtSeq->count[curSym]
    + BWTSeqTransformedOcc(bwtSeq, curSym, limits->lower);
  return limits;
}

static inline const EISeq *
BWTSeqGetEncIdxSeq(const BWTSeq *bwtSeq)
{
  assert(bwtSeq);
  return bwtSeq->seqIdx;
}

static inline struct MatchData *
EMIGetNextMatch(struct BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq)
{
  if (iter->nextMatchBWTPos < iter->bounds.lower)
  {
    iter->nextMatch.sfxArrayValue =
      BWTSeqLocateMatch(bwtSeq, iter->nextMatchBWTPos, &iter->extBits);
    {
      /* FIXME: map position back to original encoded sequence */
      iter->nextMatch.dbFile = 0;
    }
    ++iter->nextMatchBWTPos;
    return &iter->nextMatch;
  }
  else
    return NULL;
}

#endif
