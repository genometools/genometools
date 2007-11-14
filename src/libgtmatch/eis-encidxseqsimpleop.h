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
#ifndef EIS_ENCIDXSEQSIMPLEOP_H
#define EIS_ENCIDXSEQSIMPLEOP_H

#include <string.h>

#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-encidxseqpriv.h"

static inline Seqpos
EISLength(const EISeq *seq)
{
  return seq->seqLen;
}

static inline const MRAEnc *
EISGetAlphabet(const EISeq *seq)
{
  assert(seq);
  return seq->alphabet;
}

static inline Symbol
EISGetSym(EISeq *seq, Seqpos pos, EISHint hint, Env *env)
{
  assert(seq && hint && env);
  return MRAEncRevMapSymbol(seq->alphabet,
                            seq->classInfo->get(seq, pos, hint, env));
}

static inline Symbol
EISGetTransformedSym(EISeq *seq, Seqpos pos, EISHint hint, Env *env)
{
  assert(seq && hint && env);
  return seq->classInfo->get(seq, pos, hint, env);
}

static inline Seqpos
EISRank(EISeq *seq, Symbol sym, Seqpos pos, union EISHint *hint,
        Env *env)
{
  Symbol mSym;
  mSym = MRAEncMapSymbol(seq->alphabet, sym);
  return seq->classInfo->rank(seq, mSym, pos, hint, env);
}

static inline void
EISRetrieveExtraBits(EISeq *seq, Seqpos pos, int flags,
                     struct extBitsRetrieval *retval, union EISHint *hint,
                     Env *env)
{
  return seq->classInfo->expose(seq, pos, flags, retval, hint, env);
}

static inline void
initExtBitsRetrieval(struct extBitsRetrieval *r, Env *env)
{
  memset(r, 0, sizeof (struct extBitsRetrieval));
}

static inline struct extBitsRetrieval *
newExtBitsRetrieval(Env *env)
{
  struct extBitsRetrieval *retval
    = env_ma_malloc(env, sizeof (struct extBitsRetrieval));
  initExtBitsRetrieval(retval, env);
  return retval;
}

static inline void
destructExtBitsRetrieval(struct extBitsRetrieval *r, Env *env)
{
  if ((r->flags & EBRF_PERSISTENT_CWBITS) && r->cwPart)
    env_ma_free(r->cwPart, env);
  if ((r->flags & EBRF_PERSISTENT_VARBITS) && r->varPart)
    env_ma_free(r->varPart, env);
}

static inline void
deleteExtBitsRetrieval(struct extBitsRetrieval *r, Env *env)
{
  destructExtBitsRetrieval(r, env);
  env_ma_free(r, env);
}

static inline Seqpos
EISSymTransformedRank(EISeq *seq, Symbol msym, Seqpos pos,
                      union EISHint *hint, Env *env)
{
  assert(msym < MRAEncGetSize(EISGetAlphabet(seq)));
  return seq->classInfo->rank(seq, msym, pos, hint, env);
}

static inline FILE *
EISSeekToHeader(const EISeq *seqIdx, uint16_t headerID,
                uint32_t *lenRet)
{
  assert(seqIdx);
  return seqIdx->classInfo->seekToHeader(seqIdx, headerID, lenRet);
}

static inline EISHint
newEISHint(EISeq *seq, Env *env)
{
  return seq->classInfo->newHint(seq, env);
}

static inline void
deleteEISHint(EISeq *seq, EISHint hint, Env *env)
{
  return seq->classInfo->deleteHint(seq, hint, env);
}

extern int
EISPrintDiagsForPos(const EISeq *seq, Seqpos pos, FILE *fp, EISHint hint,
                    Env *env)
{
  if(seq->classInfo->printPosDiags)
    return seq->classInfo->printPosDiags(seq, pos, fp, hint, env);
  else
    return 0;
}

#endif
