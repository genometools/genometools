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

#include "libgtmatch/eis-sa-common.h"
#include "libgtmatch/eis-sa-common-priv.h"

static inline Uchar
sfxIdx2BWTSym(Seqpos sufIdx, const Encodedsequence *encseq, Readmode readmode)
{
  return sufIdx != 0
    ? getencodedchar(encseq, sufIdx - 1, readmode)
    : (Uchar) UNDEFBWTCHAR;
}

static inline size_t
EncSeqGetSubSeq(const Encodedsequence *encseq, Readmode readmode, Seqpos pos,
                size_t len, Uchar *subSeq)
{
  size_t i;
  assert(encseq);
  for (i = 0; i < len; ++i)
    subSeq[i] = getencodedchar(encseq, pos + i, readmode);
  return len;
}

static inline SeqDataReader
SASSCreateReader(SASeqSrc *src, enum sfxDataRequest request)
{
  return src->createReader(src, request);
}

static inline DefinedSeqpos
SASSGetRot0Pos(const SASeqSrc *src)
{
  return src->getRot0Pos(src);
}

static inline Seqpos
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
SASSAccessSequence(const SASeqSrc *src, Symbol *dest, Seqpos pos, size_t len)
{
  return accessSequence(src->origSequenceAccess, dest, pos, len);
}

static inline RandomSeqAccessor
SASSGetOrigSeqAccessor(const SASeqSrc *src)
{
  return src->origSequenceAccess;
}

static inline void
initSASeqSrc(SASeqSrc *src, Seqpos seqLen,
             createTranslatorFunc createTranslator,
             createReaderFunc createReader, getRot0PosFunc getRot0Pos,
             getSeqStatsFunc getSeqStats, RandomSeqAccessor origSeqAccess,
             deleteSASeqSrcFunc deleteSASS,
             newMRAEncFunc newMRAEnc,
             generatorFunc generator, void *generatorState)
{
  assert(src);
  assert(createReader || createTranslator);
  assert(getRot0Pos);
  src->seqLen = seqLen;
  src->createTranslator = createTranslator;
  if (createTranslator && !createReader)
  {
    src->createReader = SASSGenericCreateReader;
  }
  else
    src->createReader = createReader;
  src->getRot0Pos = getRot0Pos;
  src->getSeqStats = getSeqStats;
  src->origSequenceAccess = origSeqAccess;
  src->deleteSASS = deleteSASS;
  src->newMRAEnc = newMRAEnc;
  src->alphabet = NULL;
  initEmptySeqReaderSet(&src->readerSet, SFX_REQUEST_NONE, sizeof (Seqpos),
                        generator, generatorState);
  initSATaggedXltorStateList(&src->xltorStates);
}

static inline void
destructSASeqSrc(SASeqSrc *src)
{
  if (src->alphabet)
    MRAEncDelete(src->alphabet);
  destructSATaggedXltorStateList(&src->xltorStates);
  destructSeqReaderSet(&src->readerSet);
}

static inline void
SASSDelete(SASeqSrc *src)
{
  src->deleteSASS(src);
}

#endif
