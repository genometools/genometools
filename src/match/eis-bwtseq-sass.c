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

#include "core/ma_api.h"
#include "core/unused_api.h"

#include "match/eis-bwtseq.h"
#include "match/eis-bwtseq-context.h"
#include "match/eis-bwtseq-sass.h"
#include "match/eis-encidxseq.h"
#include "match/eis-list-do.h"
#include "match/eis-mrangealphabet.h"
#include "match/eis-sa-common.h"
#include "match/eis-sa-common-priv.h"

typedef struct BWTSeqReaderState BWTSeqReaderState;

struct BWTSeqReaderState
{
  BWTSeqReaderState *next;
  struct BWTSASeqSrc *backLink;
  GtUword nextReadPos;
};

struct BWTSASeqSrc
{
  struct SASeqSrc baseClass;
  const BWTSeq *bwtSeq;
  const BWTSeqContextRetriever *ctxMap;
  BWTSeqReaderState *readerStateList;
};

typedef struct BWTSASeqSrc BWTSASeqSrc;

static inline BWTSASeqSrc*
SASS2BWTSASS(SASeqSrc *baseClass)
{
  return (BWTSASeqSrc *)
    ((char *)baseClass - offsetof(BWTSASeqSrc, baseClass));
}

static inline const BWTSASeqSrc *
constSASS2BWTSASS(const SASeqSrc *baseClass)
{
  return (const BWTSASeqSrc *)
    ((const char *)baseClass - offsetof(BWTSASeqSrc, baseClass));
}

static Definedunsignedlong
BWTSASSGetRot0Pos(const struct SASeqSrc *baseClass)
{
  Definedunsignedlong val
    = { BWTSeqTerminatorPos(constSASS2BWTSASS(baseClass)->bwtSeq), true };
  return val;
}

static void
gt_deleteBWTSeqSASS(struct SASeqSrc *baseClass)
{
  destructSASeqSrc(baseClass);
  ListDo(BWTSeqReaderState, SASS2BWTSASS(baseClass)->readerStateList,
         gt_free(p));
  gt_free(SASS2BWTSASS(baseClass));
}

static size_t
BWTSASSAccessOrigSeq(const void *state,
                     Symbol *dest,
                     GtUword pos,
                     size_t len)
{
  const BWTSeqContextRetriever *ctxMap = state;
  gt_assert(state);
  gt_BWTSeqCRAccessSubseq(ctxMap, pos, len, dest);
  return len;
}

static MRAEnc *
BWTSASSNewMRAEnc(const SASeqSrc *src)
{
  const BWTSASeqSrc *bwtSASeqSrc;
  gt_assert(src);
  bwtSASeqSrc = constSASS2BWTSASS(src);
  return gt_MRAEncCopy(EISGetAlphabet(BWTSeqGetEncIdxSeq(bwtSASeqSrc->bwtSeq)));
}

static SeqDataReader
BWTSASSCreateReader(SASeqSrc *src, enum sfxDataRequest request);

SASeqSrc *
gt_BWTSeqNewSASeqSrc(const BWTSeq *bwtSeq, const BWTSeqContextRetriever *ctxMap)
{
  struct BWTSASeqSrc *newBWTSASeqSrc;
  gt_assert(bwtSeq);
  newBWTSASeqSrc = gt_malloc(sizeof (*newBWTSASeqSrc));
  {
    RandomSeqAccessor origSeqAccess;
    if (ctxMap)
    {
      origSeqAccess.accessFunc = BWTSASSAccessOrigSeq;
      origSeqAccess.state = (void *)ctxMap;
    }
    else
    {
      origSeqAccess.accessFunc = NULL;
      origSeqAccess.state = NULL;
    }
    initSASeqSrc(&newBWTSASeqSrc->baseClass,
                 BWTSeqLength(bwtSeq),
                 NULL,
                 BWTSASSCreateReader,
                 BWTSASSGetRot0Pos,
                 NULL,
                 origSeqAccess,
                 gt_deleteBWTSeqSASS,
                 BWTSASSNewMRAEnc,
                 NULL, NULL);   /* since BWTSeq can regenerate
                                 * arbitrary portions of all
                                 * associated data (when reversibly
                                 * sorted) no generator is necessary
                                 * and all readers just keep their state
                                 */
  }
  newBWTSASeqSrc->ctxMap = ctxMap;
  newBWTSASeqSrc->bwtSeq = bwtSeq;
  newBWTSASeqSrc->readerStateList = NULL;
  return &newBWTSASeqSrc->baseClass;
}

static size_t
BWTSASSReadSufTab(SeqDataSrc src, void *dest, size_t len)
{
  const BWTSeq *bwtSeq;
  struct extBitsRetrieval extBits;
  gt_assert(src);
  initExtBitsRetrieval(&extBits);
  bwtSeq = ((BWTSeqReaderState *)src)->backLink->bwtSeq;
  size_t i;
  GtUword pos = ((BWTSeqReaderState *)src)->nextReadPos;
  for (i = 0; i < len; ++i)
    ((GtUword *)dest)[i] = gt_BWTSeqLocateMatch(bwtSeq,
                                                      pos + i,
                                                      &extBits);
  ((BWTSeqReaderState *)src)->nextReadPos = pos + len;
  destructExtBitsRetrieval(&extBits);
  return len;
}

static SeqDataReader
BWTSASSMakeSufTabReader(GT_UNUSED BWTSASeqSrc *bwtSASeqSrc,
                        BWTSeqReaderState *state)
{
  SeqDataReader reader = { state, BWTSASSReadSufTab };
  return reader;
}

static size_t
BWTSASSReadBWT(SeqDataSrc src, void *dest, size_t len)
{
  const BWTSeq *bwtSeq;
  gt_assert(src);
  bwtSeq = ((BWTSeqReaderState *)src)->backLink->bwtSeq;
  size_t i;
  GtUword pos = ((BWTSeqReaderState *)src)->nextReadPos;
  for (i = 0; i < len; ++i)
    ((Symbol *)dest)[i] = BWTSeqGetSym(bwtSeq, pos + i);
  ((BWTSeqReaderState *)src)->nextReadPos = pos + len;
  return len;
}

static SeqDataReader
BWTSASSMakeBWTReader(GT_UNUSED BWTSASeqSrc *bwtSASeqSrc,
                     BWTSeqReaderState *state)
{
  SeqDataReader reader = { state, BWTSASSReadBWT };
  return reader;
}

static BWTSeqReaderState *
BWTSeqSASSAddReaderState(BWTSASeqSrc *bwtSASeqSrc)
{
  BWTSeqReaderState *newReader = gt_malloc(sizeof (*newReader));
  gt_assert(bwtSASeqSrc);
  newReader->backLink = bwtSASeqSrc;
  newReader->nextReadPos = 0;
  newReader->next = bwtSASeqSrc->readerStateList;
  bwtSASeqSrc->readerStateList = newReader;
  return newReader;
}

static SeqDataReader
BWTSASSCreateReader(SASeqSrc *src, enum sfxDataRequest rtype)
{
  struct seqDataReader reader = { NULL, NULL};
  BWTSASeqSrc *bwtSASeqSrc;
  gt_assert(src);
  bwtSASeqSrc = SASS2BWTSASS(src);
  switch (rtype)
  {
  case SFX_REQUEST_SUFTAB:
    reader = BWTSASSMakeSufTabReader(bwtSASeqSrc,
                                     BWTSeqSASSAddReaderState(bwtSASeqSrc));
    break;
  case SFX_REQUEST_BWTTAB:
    reader = BWTSASSMakeBWTReader(bwtSASeqSrc,
                                  BWTSeqSASSAddReaderState(bwtSASeqSrc));
    break;
  default:
    fprintf(stderr, "error: unimplemented request: %d, %s: %d!\n", rtype,
            __FILE__, __LINE__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  return reader;
}
