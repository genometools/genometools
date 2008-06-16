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
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"

#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseq-context.h"
#include "libgtmatch/eis-bwtseq-sass.h"
#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-list-do.h"
#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-sa-common.h"
#include "libgtmatch/eis-sa-common-priv.h"

typedef struct BWTSeqReaderState BWTSeqReaderState;

struct BWTSeqReaderState
{
  BWTSeqReaderState *next;
  struct BWTSASeqSrc *backLink;
  Seqpos nextReadPos;
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

static inline SASeqSrc *
BWTSASS2SASS(BWTSASeqSrc *bwtSASS)
{
  return &bwtSASS->baseClass;
}

static inline const SASeqSrc *
constBWTSASS2SASS(const BWTSASeqSrc *bwtSASS)
{
  return &bwtSASS->baseClass;
}

static DefinedSeqpos
BWTSASSGetRot0Pos(const struct SASeqSrc *baseClass)
{
  DefinedSeqpos val
    = { BWTSeqTerminatorPos(constSASS2BWTSASS(baseClass)->bwtSeq), true };
  return val;
}

static void
deleteBWTSeqSASS(struct SASeqSrc *baseClass)
{
  destructSASeqSrc(baseClass);
  ListDo(BWTSeqReaderState, SASS2BWTSASS(baseClass)->readerStateList,
         ma_free(p));
  ma_free(SASS2BWTSASS(baseClass));
}

static size_t
BWTSASSAccessOrigSeq(const void *state, Symbol *dest, Seqpos pos, size_t len)
{
  const BWTSeqContextRetriever *ctxMap = state;
  assert(state);
  BWTSeqCRAccessSubseq(ctxMap, pos, len, dest);
  return len;
}

static MRAEnc *
BWTSASSNewMRAEnc(const SASeqSrc *src)
{
  const BWTSASeqSrc *bwtSASeqSrc;
  assert(src);
  bwtSASeqSrc = constSASS2BWTSASS(src);
  return MRAEncCopy(EISGetAlphabet(BWTSeqGetEncIdxSeq(bwtSASeqSrc->bwtSeq)));
}

static SeqDataReader
BWTSASSCreateReader(SASeqSrc *src, enum sfxDataRequest request);

extern SASeqSrc *
BWTSeqNewSASeqSrc(const BWTSeq *bwtSeq, const BWTSeqContextRetriever *ctxMap)
{
  struct BWTSASeqSrc *newBWTSASeqSrc;
  assert(bwtSeq);
  newBWTSASeqSrc = ma_malloc(sizeof (*newBWTSASeqSrc));
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
                 deleteBWTSeqSASS,
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
BWTSASSReadSufTab(SeqDataSrc src, void *dest, size_t len, UNUSED Error *err)
{
  const BWTSeq *bwtSeq;
  struct extBitsRetrieval extBits;
  assert(src);
  initExtBitsRetrieval(&extBits);
  bwtSeq = ((BWTSeqReaderState *)src)->backLink->bwtSeq;
  size_t i;
  Seqpos pos = ((BWTSeqReaderState *)src)->nextReadPos;
  for (i = 0; i < len; ++i)
    ((Seqpos *)dest)[i] = BWTSeqLocateMatch(bwtSeq, pos + i, &extBits);
  ((BWTSeqReaderState *)src)->nextReadPos = pos + len;
  destructExtBitsRetrieval(&extBits);
  return len;
}

static SeqDataReader
BWTSASSMakeSufTabReader(UNUSED BWTSASeqSrc *bwtSASeqSrc,
                        BWTSeqReaderState *state)
{
  SeqDataReader reader = { state, BWTSASSReadSufTab };
  return reader;
}

static size_t
BWTSASSReadBWT(SeqDataSrc src, void *dest, size_t len, UNUSED Error *err)
{
  const BWTSeq *bwtSeq;
  assert(src);
  bwtSeq = ((BWTSeqReaderState *)src)->backLink->bwtSeq;
  size_t i;
  Seqpos pos = ((BWTSeqReaderState *)src)->nextReadPos;
  for (i = 0; i < len; ++i)
    ((Symbol *)dest)[i] = BWTSeqGetSym(bwtSeq, pos + i);
  ((BWTSeqReaderState *)src)->nextReadPos = pos + len;
  return len;
}

static SeqDataReader
BWTSASSMakeBWTReader(UNUSED BWTSASeqSrc *bwtSASeqSrc, BWTSeqReaderState *state)
{
  SeqDataReader reader = { state, BWTSASSReadBWT };
  return reader;
}

static BWTSeqReaderState *
BWTSeqSASSAddReaderState(BWTSASeqSrc *bwtSASeqSrc)
{
  BWTSeqReaderState *newReader = ma_malloc(sizeof (*newReader));
  assert(bwtSASeqSrc);
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
  assert(src);
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
  case SFX_REQUEST_LCPTAB:
  default:
    fprintf(stderr, "error: unimplemented request: %d, %s: %d!\n", rtype,
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  return reader;
}
