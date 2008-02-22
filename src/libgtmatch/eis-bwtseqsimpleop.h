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

static inline struct SeqposPair
BWTSeqTransformedPosPairOcc(const BWTSeq *bwtSeq, Symbol tSym,
                            Seqpos posA, Seqpos posB)
{
  assert(bwtSeq);
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
/*       - ((pos > bwtSeq->longest)?1:0); */
  else /* tSym == not flattened terminator == alphabetSize - 1 */
  {
    struct SeqposPair occ;
    assert(tSym == bwtSeq->alphabetSize - 1);
    occ.a = (posA > bwtSeq->longest)?1:0;
    occ.b = (posB > bwtSeq->longest)?1:0;
    return occ;
  }
}

static inline Seqpos
BWTSeqOcc(const BWTSeq *bwtSeq, Symbol sym, Seqpos pos)
{
  assert(bwtSeq);
  Symbol tSym = MRAEncMapSymbol(bwtSeq->alphabet, sym);
  return BWTSeqTransformedOcc(bwtSeq, tSym, pos);
}

static inline struct SeqposPair
BWTSeqPosPairOcc(const BWTSeq *bwtSeq, Symbol sym, Seqpos posA, Seqpos posB)
{
  assert(bwtSeq);
  Symbol tSym = MRAEncMapSymbol(bwtSeq->alphabet, sym);
  return BWTSeqTransformedPosPairOcc(bwtSeq, tSym, posA, posB);
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
  assert(limits->start < bwtSeq->count[bwtSeq->alphabetSize]
         && limits->end < bwtSeq->count[bwtSeq->alphabetSize]);
  alphabet = BWTSeqGetAlphabet(bwtSeq);
  curSym = MRAEncMapSymbol(alphabet, nextSym);
  assert(MRAEncSymbolHasValidMapping(alphabet, curSym));
  limits->start = bwtSeq->count[curSym]
    + BWTSeqTransformedOcc(bwtSeq, curSym, limits->start);
  limits->end = bwtSeq->count[curSym]
    + BWTSeqTransformedOcc(bwtSeq, curSym, limits->end);
  return limits;
}

static inline const EISeq *
BWTSeqGetEncIdxSeq(const BWTSeq *bwtSeq)
{
  assert(bwtSeq);
  return bwtSeq->seqIdx;
}

static inline bool
EMIGetNextMatch(struct BWTSeqExactMatchesIterator *iter, Seqpos *pos,
                const BWTSeq *bwtSeq)
{
  if (iter->nextMatchBWTPos < iter->bounds.end)
  {
    *pos = BWTSeqLocateMatch(bwtSeq, iter->nextMatchBWTPos, &iter->extBits);
    ++iter->nextMatchBWTPos;
    return true;
  }
  else
    return false;
}

static inline Seqpos pckfindfirstmatch(const BWTSeq *bwtSeq,Seqpos lowerbound)
{
  struct extBitsRetrieval extBits;
  Seqpos pos;

  initExtBitsRetrieval(&extBits);
  pos = BWTSeqLocateMatch(bwtSeq,lowerbound,&extBits);
  destructExtBitsRetrieval(&extBits);
  return pos;
}

#endif
