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

#include <stdlib.h>
#include <string.h>

#include "core/ma_api.h"
#include "core/unused_api.h"
#include "match/eis-sequencemultiread.h"
#include "match/sarr-def.h"
#include "match/eis-list-do.h"

struct seqReaderState
{
  struct seqReaderState *next;
  GtUword nextReadPos;           /**< next position to be read/consumed */
  int tag;                      /**< set of bits defined by
                                 * application, can be used for value
                                 * or as bit set */
  SeqReaderSet *readerSet;
  consumerID id;
  void *dest;
  SeqDataTranslator xltor;
};

struct seqSinkState
{
  int tag;                      /**< set of bits defined by
                                 * application, can be used for value
                                 * or as bit set */
  consumerID id;
  SeqDataWriter writer;
};

static inline GtUword
seqReaderSetGetConsumerNextPos(struct seqReaderState *consumer)
{
  gt_assert(consumer);
  return consumer->nextReadPos;
}

static inline void
seqReaderSetSetConsumerNextPos(struct seqReaderState *consumer,
                               GtUword newPos)
{
  gt_assert(consumer && consumer->nextReadPos <= newPos);
  consumer->nextReadPos = newPos;
}

static void
seqReaderSetMove2Backlog(void *backlogState, const void *seqData,
                         GtUword requestStart, size_t requestLen);

static size_t
seqReaderSetRead(void *state, void *dest, size_t len);

/*
int
gt_initSeqReaderSet(SeqReaderSet *readerSet, int initialSuperSet,
                 int numConsumers, int *tags, SeqDataTranslator xltors[],
                 SeqDataReader *generatedReaders, size_t seqElemSize,
                 generatorFunc generator, void *generatorState)
{
  int i;
  gt_initEmptySeqReaderSet(readerSet, initialSuperSet, seqElemSize,
                        generator, generatorState);
  for (i = 0; i < numConsumers; ++i)
  {
    generatedReaders[i]
      = gt_seqReaderSetRegisterConsumer(readerSet, tags[i], xltors[i]);
    if (!generatedReaders[i].readData)
      break;
  }
  return i;
}
*/

void
gt_initEmptySeqReaderSet(SeqReaderSet *readerSet, int initialSuperSet,
                         bool fromSuffixsortspace,
                         size_t seqElemSize, generatorFunc generator,
                         void *generatorState)
{
  readerSet->fromSuffixsortspace = fromSuffixsortspace;
  readerSet->tagSuperSet = initialSuperSet;
  readerSet->numAutoConsumers = readerSet->numConsumers = 0;
  readerSet->consumerList = NULL;
  readerSet->autoConsumerList = NULL;
  readerSet->backlogStartPos = 0;
  readerSet->backlogLen = readerSet->backlogSize = 0;
  readerSet->backlogElemSize = seqElemSize;
  readerSet->seqDataBacklog = NULL; /* only allocate backlog when it's
                                     * really needed! */
  readerSet->generatorState = generatorState;
  readerSet->generator = generator;
}

void
gt_destructSeqReaderSet(SeqReaderSet *readerSet)
{
  gt_assert(readerSet);
  ListDo(struct seqReaderState, readerSet->consumerList, gt_free(p));
  if (readerSet->autoConsumerList)
    gt_free(readerSet->autoConsumerList);
  if (readerSet->seqDataBacklog)
    gt_free(readerSet->seqDataBacklog);
}

SeqDataReader
gt_seqReaderSetRegisterConsumer(SeqReaderSet *readerSet, int tag,
                             SeqDataTranslator xltor)
{
  int availId = readerSet->numConsumers++;
  SeqDataReader reader;
  struct seqReaderState *state;
  gt_assert(readerSet);
  state = gt_malloc(sizeof (*state));
  state->readerSet = readerSet;
  state->id = availId;
  state->xltor = xltor;
  state->nextReadPos = 0;
  state->tag = tag;
  state->next = readerSet->consumerList;
  readerSet->consumerList = state;

  readerSet->tagSuperSet |= tag;

  reader.src = state;
  reader.readData = seqReaderSetRead;
  return reader;
}

bool
gt_seqReaderSetRegisterAutoConsumer(SeqReaderSet *readerSet, int tag,
                                 SeqDataWriter writer)
{
  int availId = readerSet->numConsumers++;
  gt_assert(false);
  gt_assert(readerSet);
  {
    struct seqSinkState *temp;
    if (!(temp = gt_realloc(readerSet->autoConsumerList,
                            sizeof (*temp) * (availId + 1))))
      return false;
    readerSet->autoConsumerList = temp;
  }
  {
    struct seqSinkState *state = readerSet->autoConsumerList + availId;
    state->id = availId;
    state->writer = writer;
    state->tag = tag;

    readerSet->tagSuperSet |= tag;
  }
  return true;
}

static size_t
seqReaderSetRead(void *src, void *dest, size_t len)
{
  struct seqReaderState *state = src;
  struct seqReaderSet *readerSet;
  GtUword pos;
  size_t elemsLeft;
  gt_assert(src);
  readerSet = state->readerSet;
  elemsLeft = len;
  pos = seqReaderSetGetConsumerNextPos(state);
  while (elemsLeft)
  {
    if (pos >= readerSet->backlogStartPos + readerSet->backlogLen)
    {
      /* pos is in block that's governed by generator */
      size_t elemsGenerated =
        readerSet->generator(readerSet->generatorState, readerSet,
                             seqReaderSetMove2Backlog, dest, pos, elemsLeft,
                             state->xltor);
      seqReaderSetSetConsumerNextPos(state, pos += elemsGenerated);
      elemsLeft -= elemsGenerated;
      break;
    }
    else if (pos >= readerSet->backlogStartPos
             /* implicit because of above condition
                && pos < readerSet->backlogStartPos + readerSet->backlogLen */
      )
    {
      /* pos is in backlog */
      GtUword subLen, charsWritten;

      subLen = MIN(elemsLeft, readerSet->backlogStartPos - pos
                              + readerSet->backlogLen);
      gt_assert(state->xltor.translateData == NULL);
      /*
      charsWritten = SDRTranslate(state->xltor, dest,
                                  (char *)readerSet->seqDataBacklog
                                  + (pos - readerSet->backlogStartPos)
                                  * readerSet->backlogElemSize,
                                  subLen);
      */
      memcpy(dest, (char *)readerSet->seqDataBacklog
                   + (pos - readerSet->backlogStartPos)
                   * readerSet->backlogElemSize,
             subLen * state->xltor.state.elemSize);
      charsWritten = subLen * state->xltor.state.elemSize;
      elemsLeft -= subLen;
      seqReaderSetSetConsumerNextPos(state, pos += subLen);
      dest = (char *)dest + charsWritten;
    }
    else
    {
      fprintf(stderr, "fatal error at file %s, line %d\n", __FILE__, __LINE__);
      abort();
    }
  }
  return len - elemsLeft;
}

#if 0
/* These might be needed later when translation results values are
 * kept in backlog too */
static inline GtUword
seqReaderSetFindMinOpenRequestBySet(SeqReaderSet *readerSet, int setSpec)
{
  GtUword min = -1;
  int i;
  gt_assert(readerSet);
  for (i = 0; i < readerSet->numConsumers; ++i)
    if (readerSet->consumers[i].tag & setSpec)
      min = MIN(min, readerSet->consumers[i].nextReadPos);
  return min;
}

static inline GtUword
seqReaderSetFindMinOpenRequestByVal(SeqReaderSet *readerSet, int tagVal)
{
  GtUword min = -1;
  int i;
  gt_assert(readerSet);
  for (i = 0; i < readerSet->numConsumers; ++i)
    if (readerSet->consumers[i].tag == tagVal)
      min = MIN(min, readerSet->consumers[i].nextReadPos);
  return min;
}
#endif

static inline GtUword
seqReaderSetFindMinOpenRequest(SeqReaderSet *readerSet)
{
  struct seqReaderState *p;
  GtUword min = -1;
  gt_assert(readerSet);
  p = readerSet->consumerList;
  while (p)
  {
    min = MIN(min, seqReaderSetGetConsumerNextPos(p));
    p = p->next;
  }
  return min;
}

static void
seqReaderSetMove2Backlog(void *backlogState, const void *seqData,
                         GtUword requestStart, size_t requestLen)
{
  GtUword requestMinPos;
  struct seqReaderSet *readerSet = backlogState;
  gt_assert(backlogState && (requestLen?(seqData!=NULL):1));
  requestMinPos = seqReaderSetFindMinOpenRequest(readerSet);
  /* 1. pass all data to be invalidated to automatic sinks */
  {
    int i, numAutoConsumers = readerSet->numAutoConsumers;
    struct seqSinkState *sinks = readerSet->autoConsumerList;
    /* The following must hold, as writer is not defined */
    gt_assert(numAutoConsumers == 0);
    for (i = 0; i < numAutoConsumers; ++i)
    {
      SDWWrite(sinks[i].writer, seqData, requestLen);
    }
  }
  /* 2. move still unread old values as far as possible to head of copy */
  gt_assert(requestMinPos >= readerSet->backlogStartPos);
  if (requestMinPos < readerSet->backlogStartPos + readerSet->backlogLen)
  {
    size_t backlogUnread
      = readerSet->backlogLen - requestMinPos + readerSet->backlogStartPos;
    memmove(readerSet->seqDataBacklog,
            (char *)readerSet->seqDataBacklog
            + (requestMinPos - readerSet->backlogStartPos)
            * readerSet->backlogElemSize,
            readerSet->backlogElemSize * backlogUnread);
    readerSet->backlogLen = backlogUnread;
    readerSet->backlogStartPos = requestMinPos;
  }
  else
  {
    readerSet->backlogLen = 0;
    readerSet->backlogStartPos = requestMinPos;
  }
  {
    GtUword copyStartPos = MAX(requestMinPos, requestStart);
    size_t copyLen = requestLen - copyStartPos + requestStart;
    /* 3. extend backlog to also accept all invalidated values still needed */
    if (copyLen)
    {
      GtUword *destSptr;

      size_t idx, backlogSizeLeft
        = readerSet->backlogSize - readerSet->backlogLen;
      if (copyLen > backlogSizeLeft)
      {
        size_t newSize = readerSet->backlogLen + copyLen;
        readerSet->seqDataBacklog
          = gt_realloc(readerSet->seqDataBacklog, readerSet->backlogElemSize
                       * newSize);
        readerSet->backlogSize = newSize;
      }
      destSptr = ((GtUword *) readerSet->seqDataBacklog) +
                 readerSet->backlogLen;
      if (readerSet->fromSuffixsortspace)
      {
        GtSuffixsortspace *srcSptr = (GtSuffixsortspace *) seqData;
        for (idx = 0; idx< copyLen; idx++)
        {
          destSptr[idx]
            = gt_suffixsortspace_getdirect(srcSptr,
                                           copyStartPos - requestStart + idx);
        }
      } else
      {
        ESASuffixptr *srcSptr
          = ((ESASuffixptr *) seqData) + (copyStartPos - requestStart);
        for (idx = 0; idx< copyLen; idx++)
        {
          destSptr[idx] = ESASUFFIXPTRGET(srcSptr,idx);
        }
      }
      /*
      memcpy((char *)readerSet->seqDataBacklog
             + readerSet->backlogLen * readerSet->backlogElemSize,
             (char *)seqData
             + (copyStartPos - requestStart) * readerSet->backlogElemSize,
             readerSet->backlogElemSize * copyLen);
      */
      readerSet->backlogLen += copyLen;
    }
  }
}