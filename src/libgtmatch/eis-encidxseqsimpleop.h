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
#include "libgtcore/ma.h"
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
EISGetSym(EISeq *seq, Seqpos pos, EISHint hint)
{
  assert(seq && hint);
  return MRAEncRevMapSymbol(seq->alphabet,
                            seq->classInfo->get(seq, pos, hint));
}

static inline Symbol
EISGetTransformedSym(EISeq *seq, Seqpos pos, EISHint hint)
{
  assert(seq && hint);
  return seq->classInfo->get(seq, pos, hint);
}

static inline Seqpos
EISRank(EISeq *seq, Symbol sym, Seqpos pos, union EISHint *hint)
{
  Symbol mSym;
  mSym = MRAEncMapSymbol(seq->alphabet, sym);
  return seq->classInfo->rank(seq, mSym, pos, hint);
}

static inline void
EISRetrieveExtraBits(EISeq *seq, Seqpos pos, int flags,
                     struct extBitsRetrieval *retval, union EISHint *hint)
{
  return seq->classInfo->expose(seq, pos, flags, retval, hint);
}

static inline void
initExtBitsRetrieval(struct extBitsRetrieval *r)
{
  memset(r, 0, sizeof (struct extBitsRetrieval));
}

static inline struct extBitsRetrieval *
newExtBitsRetrieval()
{
  struct extBitsRetrieval *retval = ma_malloc(sizeof (struct extBitsRetrieval));
  initExtBitsRetrieval(retval);
  return retval;
}

static inline void
destructExtBitsRetrieval(struct extBitsRetrieval *r)
{
  if ((r->flags & EBRF_PERSISTENT_CWBITS) && r->cwPart)
    ma_free(r->cwPart);
  if ((r->flags & EBRF_PERSISTENT_VARBITS) && r->varPart)
    ma_free(r->varPart);
}

static inline void
deleteExtBitsRetrieval(struct extBitsRetrieval *r)
{
  destructExtBitsRetrieval(r);
  ma_free(r);
}

static inline Seqpos
EISSymTransformedRank(EISeq *seq, Symbol msym, Seqpos pos,
                      union EISHint *hint)
{
  assert(msym < MRAEncGetSize(EISGetAlphabet(seq)));
  return seq->classInfo->rank(seq, msym, pos, hint);
}

static inline FILE *
EISSeekToHeader(const EISeq *seqIdx, uint16_t headerID,
                uint32_t *lenRet)
{
  assert(seqIdx);
  return seqIdx->classInfo->seekToHeader(seqIdx, headerID, lenRet);
}

static inline EISHint
newEISHint(const EISeq *seq)
{
  return seq->classInfo->newHint(seq);
}

static inline void
deleteEISHint(EISeq *seq, EISHint hint)
{
  return seq->classInfo->deleteHint(seq, hint);
}

static inline int
EISPrintDiagsForPos(const EISeq *seq, Seqpos pos, FILE *fp, EISHint hint)
{
  if (seq->classInfo->printPosDiags)
    return seq->classInfo->printPosDiags(seq, pos, fp, hint);
  else
    return 0;
}

static inline int
EISPrintExtDiagsForPos(const EISeq *seq, Seqpos pos, FILE *fp, EISHint hint)
{
  if (seq->classInfo->printExtPosDiags)
    return seq->classInfo->printExtPosDiags(seq, pos, fp, hint);
  else
    return 0;
}

#endif
