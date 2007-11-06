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
  return EISGetAlphabet(seq->seqIdx);
}

static inline Seqpos
BWTSeqLength(const BWTSeq *seq)
{
  return EISLength(seq->seqIdx);
}

static inline Seqpos
BWTSeqTransformedOcc(const BWTSeq *bwtSeq, Symbol tsym, Seqpos pos, Env *env)
{
  assert(bwtSeq && env);
  return EISSymTransformedRank(bwtSeq->seqIdx, tsym, pos, bwtSeq->hint, env);
}

static inline Seqpos
BWTSeqOcc(const BWTSeq *bwtSeq, Symbol sym, Seqpos pos, Env *env)
{
  assert(bwtSeq && env);
  return EISRank(bwtSeq->seqIdx, sym, pos, bwtSeq->hint, env);
}

static inline Seqpos
BWTSeqLFMap(const BWTSeq *bwtSeq, Seqpos pos, Env *env)
{
  Symbol tSym = EISGetTransformedSym(bwtSeq->seqIdx, pos, bwtSeq->hint, env);
  return bwtSeq->count[tSym] + BWTSeqTransformedOcc(bwtSeq, tSym, pos, env);
}

static inline Seqpos
BWTSeqAggCount(const BWTSeq *bwtSeq, Symbol sym, Env *env)
{
  Symbol  tSym;
  assert(bwtSeq && env);
  assert(MRAEncSymbolHasValidMapping(BWTSeqGetAlphabet(bwtSeq), sym));
  tSym = MRAEncMapSymbol(BWTSeqGetAlphabet(bwtSeq), sym);
  return bwtSeq->count[tSym];
}

static inline Seqpos
BWTSeqAggTransformedCount(const BWTSeq *bwtSeq, Symbol tSym, Env *env)
{
  assert(bwtSeq && env);
  assert(tSym <= bwtSeq->alphabetSize);
  return bwtSeq->count[tSym];
}

static struct matchBound *
BWTSeqIncrMatch(const BWTSeq *bwtSeq, struct matchBound *limits,
                Symbol nextSym, Env *env)
{
  const MRAEnc *alphabet;
  Symbol curSym;
  assert(bwtSeq && limits && env);
  assert(limits->upper < bwtSeq->count[bwtSeq->alphabetSize]
         && limits->lower < bwtSeq->count[bwtSeq->alphabetSize]);
  alphabet = EISGetAlphabet(bwtSeq->seqIdx);
  curSym = MRAEncMapSymbol(alphabet, nextSym);
  assert(MRAEncSymbolHasValidMapping(alphabet, curSym));
  limits->upper = bwtSeq->count[curSym]
    + BWTSeqTransformedOcc(bwtSeq, curSym, limits->upper, env);
  limits->lower = bwtSeq->count[curSym]
    + BWTSeqTransformedOcc(bwtSeq, curSym, limits->lower, env);
  return limits;
}

static inline const EISeq *
BWTSeqGetEncIdxSeq(const BWTSeq *bwtSeq)
{
  assert(bwtSeq);
  return bwtSeq->seqIdx;
}

#endif
