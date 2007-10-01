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

#include "mrangealphabet.h"
#include "encidxseq.h"
#include "encidxseqpriv.h"

staticifinline inline Seqpos
EISLength(struct encIdxSeq *seq)
{
  return seq->seqLen;
}

staticifinline inline const MRAEnc *
EISGetAlphabet(const struct encIdxSeq *seq)
{
  assert(seq);
  return seq->alphabet;
}

staticifinline inline Symbol
EISGetSym(struct encIdxSeq *seq, Seqpos pos, union EISHint *hint, Env *env)
{
  assert(seq && hint && env);
  return MRAEncRevMapSymbol(seq->alphabet,
                            seq->classInfo->get(seq, pos, hint, env));
}

staticifinline inline Symbol
EISGetTransformedSym(struct encIdxSeq *seq, Seqpos pos, EISHint hint, Env *env)
{
  assert(seq && hint && env);
  return seq->classInfo->get(seq, pos, hint, env);
}

staticifinline inline Seqpos
EISRank(struct encIdxSeq *seq, Symbol sym, Seqpos pos, union EISHint *hint,
        Env *env)
{
  Symbol mSym;
  mSym = MRAEncMapSymbol(seq->alphabet, sym);
  return seq->classInfo->rank(seq, mSym, pos, hint, env);
}

staticifinline inline void
EISRetrieveExtraBits(struct encIdxSeq *seq, Seqpos pos, int flags,
                     struct extBitsRetrieval *retval, union EISHint *hint,
                     Env *env)
{
  return seq->classInfo->expose(seq, pos, flags, retval, hint, env);
}

staticifinline inline void
initExtBitsRetrieval(struct extBitsRetrieval *r, Env *env)
{
  memset(r, 0, sizeof(struct extBitsRetrieval));
}

staticifinline inline struct extBitsRetrieval *
newExtBitsRetrieval(Env *env)
{
  struct extBitsRetrieval *retval 
    = env_ma_malloc(env, sizeof(struct extBitsRetrieval));
  initExtBitsRetrieval(retval, env);
  return retval;
}

staticifinline inline void
destructExtBitsRetrieval(struct extBitsRetrieval *r, Env *env)
{
  if((r->flags & EBRF_PERSISTENT_CWBITS) && r->cwPart)
    env_ma_free(r->cwPart, env);
  if((r->flags & EBRF_PERSISTENT_VARBITS) && r->varPart)
    env_ma_free(r->varPart, env);
}

staticifinline inline void
deleteExtBitsRetrieval(struct extBitsRetrieval *r, Env *env)
{
  destructExtBitsRetrieval(r, env);
  env_ma_free(r, env);
}

staticifinline inline Seqpos
EISSymTransformedRank(struct encIdxSeq *seq, Symbol msym, Seqpos pos,
                      union EISHint *hint, Env *env)
{
  assert(msym < MRAEncGetSize(EISGetAlphabet(seq)));
  return seq->classInfo->rank(seq, msym, pos, hint, env);
}

staticifinline inline FILE *
EISSeekToHeader(const struct encIdxSeq *seqIdx, uint16_t headerID,
                uint32_t *lenRet)
{
  assert(seqIdx);
  return seqIdx->classInfo->seekToHeader(seqIdx, headerID, lenRet);
}

staticifinline inline EISHint
newEISHint(struct encIdxSeq *seq, Env *env)
{
  return seq->classInfo->newHint(seq, env);
}

staticifinline inline void
deleteEISHint(struct encIdxSeq *seq, EISHint hint, Env *env)
{
  return seq->classInfo->deleteHint(seq, hint, env);
}


