/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/

#include "bwtseq.h"
#include "bwtseqpriv.h"

staticifinline inline const MRAEnc *
BWTSeqGetAlphabet(const BWTSeq *seq)
{
  return EISGetAlphabet(seq->seqIdx);
}

staticifinline inline Seqpos
BWTSeqLength(const BWTSeq *seq)
{
  return EISLength(seq->seqIdx);
}

staticifinline inline Seqpos
BWTSeqTransformedOcc(const BWTSeq *bwtSeq, Symbol tsym, Seqpos pos, Env *env)
{
  assert(bwtSeq && env);
  return EISSymTransformedRank(bwtSeq->seqIdx, tsym, pos, bwtSeq->hint, env);
}

staticifinline inline Seqpos
BWTSeqOcc(const BWTSeq *bwtSeq, Symbol sym, Seqpos pos, Env *env)
{
  assert(bwtSeq && env);
  return EISRank(bwtSeq->seqIdx, sym, pos, bwtSeq->hint, env);
}

staticifinline inline Seqpos
BWTSeqLFMap(const BWTSeq *bwtSeq, Seqpos pos, Env *env)
{
  Symbol tSym = EISGetTransformedSym(bwtSeq->seqIdx, pos, bwtSeq->hint, env);
  return bwtSeq->count[tSym] + BWTSeqTransformedOcc(bwtSeq, tSym, pos, env);
}

staticifinline inline Seqpos
BWTSeqAggCount(const BWTSeq *bwtSeq, Symbol sym, Env *env)
{
  Symbol  tSym;
  assert(bwtSeq && env);
  assert(MRAEncSymbolHasValidMapping(BWTSeqGetAlphabet(bwtSeq), sym));
  tSym = MRAEncMapSymbol(BWTSeqGetAlphabet(bwtSeq), sym);
  return bwtSeq->count[tSym];
}

staticifinline inline Seqpos
BWTSeqAggTransformedCount(const BWTSeq *bwtSeq, Symbol tSym, Env *env)
{
  assert(bwtSeq && env);
  assert(tSym <= bwtSeq->alphabetSize);
  return bwtSeq->count[tSym];
}

staticifinline struct matchBound *
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

staticifinline inline const EISeq *
BWTSeqGetEncIdxSeq(const BWTSeq *bwtSeq)
{
  assert(bwtSeq);
  return bwtSeq->seqIdx;
}


