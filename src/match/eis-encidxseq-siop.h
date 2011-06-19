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

#ifndef EIS_ENCIDXSEQ_SIOP_H
#define EIS_ENCIDXSEQ_SIOP_H

#include <string.h>
#include "core/ma_api.h"
#include "match/eis-mrangealphabet.h"
#include "match/eis-encidxseq.h"
#include "match/eis-encidxseq-priv.h"

static inline unsigned long
EISLength(const EISeq *seq)
{
  return seq->seqLen;
}

static inline const MRAEnc *
EISGetAlphabet(const EISeq *seq)
{
  gt_assert(seq);
  return seq->alphabet;
}

static inline Symbol
EISGetSym(EISeq *seq, unsigned long pos, EISHint hint)
{
  gt_assert(seq && hint);
  return MRAEncRevMapSymbol(seq->alphabet,
                            seq->classInfo->get(seq, pos, hint));
}

static inline Symbol
EISGetTransformedSym(EISeq *seq, unsigned long pos, EISHint hint)
{
  gt_assert(seq && hint);
  return seq->classInfo->get(seq, pos, hint);
}

static inline unsigned long
EISRank(EISeq *seq, Symbol sym, unsigned long pos, union EISHint *hint)
{
  Symbol mSym;
  mSym = MRAEncMapSymbol(seq->alphabet, sym);
  return seq->classInfo->rank(seq, mSym, pos, hint);
}

static inline unsigned long
EISSymTransformedRank(EISeq *seq, Symbol tSym, unsigned long pos,
                      union EISHint *hint)
{
  gt_assert(tSym < gt_MRAEncGetSize(EISGetAlphabet(seq)));
  return seq->classInfo->rank(seq, tSym, pos, hint);
}

static inline GtUlongPair
EISPosPairRank(EISeq *seq, Symbol sym, unsigned long posA, unsigned long posB,
               union EISHint *hint)
{
  Symbol tSym;
  tSym = MRAEncMapSymbol(seq->alphabet, sym);
  return seq->classInfo->posPairRank(seq, tSym, posA, posB, hint);
}

static inline void
EISRangeRank(EISeq *seq, AlphabetRangeID range, unsigned long pos,
             unsigned long *rankCounts, union EISHint *hint)
{
  return seq->classInfo->rangeRank(seq, range, pos, rankCounts, hint);
}

static inline void
EISPosPairRangeRank(EISeq *seq, AlphabetRangeID range, unsigned long posA,
                    unsigned long posB, unsigned long *rankCounts,
                    union EISHint *hint)
{
  seq->classInfo->posPairRangeRank(seq, range, posA, posB, rankCounts, hint);
}

static inline GtUlongPair
EISSymTransformedPosPairRank(EISeq *seq, Symbol tSym, unsigned long posA,
                             unsigned long posB, union EISHint *hint)
{
  if (tSym >= gt_MRAEncGetSize(EISGetAlphabet(seq)))
  {
    fprintf(stderr,"tsym=%lu,gt_MRAEncGetSize(EISGetAlphabet(seq)=%lu\n",
           (unsigned long) tSym,
           (unsigned long) gt_MRAEncGetSize(EISGetAlphabet(seq)));
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_assert(tSym < gt_MRAEncGetSize(EISGetAlphabet(seq)));
  return seq->classInfo->posPairRank(seq, tSym, posA, posB, hint);
}

static inline void
EISRetrieveExtraBits(EISeq *seq, unsigned long pos, int flags,
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
  struct extBitsRetrieval *retval
    = gt_malloc(sizeof (struct extBitsRetrieval));
  initExtBitsRetrieval(retval);
  return retval;
}

static inline void
destructExtBitsRetrieval(struct extBitsRetrieval *r)
{
  if ((r->flags & EBRF_PERSISTENT_CWBITS) && r->cwPart)
    gt_free(r->cwPart);
  if ((r->flags & EBRF_PERSISTENT_VARBITS) && r->varPart)
    gt_free(r->varPart);
}

static inline void
deleteExtBitsRetrieval(struct extBitsRetrieval *r)
{
  destructExtBitsRetrieval(r);
  gt_free(r);
}

static inline FILE *
EISSeekToHeader(const EISeq *seqIdx, uint16_t headerID,
                uint32_t *lenRet)
{
  gt_assert(seqIdx);
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
EISPrintDiagsForPos(const EISeq *seq, unsigned long pos, FILE *fp, EISHint hint)
{
  if (seq->classInfo->printPosDiags)
    return seq->classInfo->printPosDiags(seq, pos, fp, hint);
  else
    return 0;
}

static inline int
EISPrintExtDiagsForPos(const EISeq *seq, unsigned long pos, FILE *fp,
                       EISHint hint)
{
  if (seq->classInfo->printExtPosDiags)
    return seq->classInfo->printExtPosDiags(seq, pos, fp, hint);
  else
    return 0;
}

#endif
