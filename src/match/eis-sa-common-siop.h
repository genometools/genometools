/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#ifndef EIS_SA_COMMON_SIOP_H
#define EIS_SA_COMMON_SIOP_H

#include "match/eis-sa-common.h"
#include "match/eis-sa-common-priv.h"

static inline GtUchar
sfxIdx2BWTSym(GtUword sufIdx, const GtEncseq *encseq,
              GtReadmode readmode)
{
  return sufIdx != 0
    ? gt_encseq_get_encoded_char(encseq, sufIdx - 1, readmode)
    : (GtUchar) UNDEFBWTCHAR;
}

static inline size_t
EncSeqGetSubSeq(const GtEncseq *encseq, GtReadmode readmode,
                GtUword pos, size_t len, GtUchar *subSeq)
{
  size_t i;
  gt_assert(encseq);
  for (i = 0; i < len; ++i)
    subSeq[i] = gt_encseq_get_encoded_char(encseq, pos + i, readmode);
  return len;
}

static inline SeqDataReader
SASSCreateReader(SASeqSrc *src, enum sfxDataRequest request)
{
  return src->createReader(src, request);
}

static inline Definedunsignedlong
SASSGetRot0Pos(const SASeqSrc *src)
{
  return src->getRot0Pos(src);
}

static inline GtUword
SASSGetLength(const SASeqSrc *src)
{
  return src->seqLen;
}

static inline MRAEnc *
SASSNewMRAEnc(const SASeqSrc *src)
{
  return src->newMRAEnc(src);
}

static inline const MRAEnc *
SASSGetMRAEnc(SASeqSrc *src)
{
  return (src->alphabet ?
          src->alphabet : (src->alphabet = SASSNewMRAEnc(src)));
}

static inline const struct seqStats *
SASSGetSeqStats(const SASeqSrc *src)
{
  if (src->getSeqStats)
    return src->getSeqStats(src);
  else
    return NULL;
}

static inline size_t
SASSAccessSequence(const SASeqSrc *src, Symbol *dest,
                   GtUword pos, size_t len)
{
  return accessSequence(src->origSequenceAccess, dest, pos, len);
}

static inline RandomSeqAccessor
SASSGetOrigSeqAccessor(const SASeqSrc *src)
{
  return src->origSequenceAccess;
}

static inline void
initSASeqSrc(SASeqSrc *src, GtUword seqLen,
             createTranslatorFunc createTranslator,
             createReaderFunc createReader, getRot0PosFunc getRot0Pos,
             getSeqStatsFunc getSeqStats, RandomSeqAccessor origSeqAccess,
             deleteSASeqSrcFunc deleteSASS,
             newMRAEncFunc newMRAEnc,
             generatorFunc generator, void *generatorState)
{
  gt_assert(src);
  gt_assert(createReader || createTranslator);
  gt_assert(getRot0Pos);
  src->seqLen = seqLen;
  src->createTranslator = createTranslator;
  /* createTranslator is not NULL iff read from suffixerator */
  /* createReader is not NULL iff read from suffixarray */
  if (createTranslator && !createReader)
  {
    src->createReader = gt_SASSGenericCreateReader;
  } else
  {
    src->createReader = createReader;
  }
  src->getRot0Pos = getRot0Pos;
  src->getSeqStats = getSeqStats;
  src->origSequenceAccess = origSeqAccess;
  src->deleteSASS = deleteSASS;
  src->newMRAEnc = newMRAEnc;
  src->alphabet = NULL;
  gt_initEmptySeqReaderSet(&src->readerSet,
                           SFX_REQUEST_NONE,
                           (createTranslator != NULL) ? true : false,
                           sizeof (GtUword),
                           generator, generatorState);
  gt_initSATaggedXltorStateList(&src->xltorStates);
}

static inline void
destructSASeqSrc(SASeqSrc *src)
{
  if (src->alphabet)
    gt_MRAEncDelete(src->alphabet);
  gt_destructSATaggedXltorStateList(&src->xltorStates);
  gt_destructSeqReaderSet(&src->readerSet);
}

static inline void
SASSDelete(SASeqSrc *src)
{
  src->deleteSASS(src);
}

#endif
