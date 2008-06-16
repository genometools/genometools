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

#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtmatch/eis-sequencemultiread.h"
#include "libgtmatch/eis-list-do.h"

struct seqReaderState
{
  struct seqReaderState *next;
  Seqpos nextReadPos;           /**< next position to be read/consumed */
  int tag;                      /**< set of bits defined by
                                 * application, can be used for value
                                 * or as bit set */
  SeqReaderSet *readerSet;
  consumerID id;
  SeqDataTranslator xltor;
};

static inline Seqpos
seqReaderSetGetConsumerNextPos(struct seqReaderState *consumer)
{
  assert(consumer);
  return consumer->nextReadPos;
}

static inline void
seqReaderSetSetConsumerNextPos(struct seqReaderState *consumer,
                               Seqpos newPos)
{
  assert(consumer && consumer->nextReadPos <= newPos);
  consumer->nextReadPos = newPos;
}

static void
seqReaderSetMove2Backlog(void *backlogState, const void *seqData,
                         Seqpos requestStart, size_t requestLen);

static size_t
seqReaderSetRead(void *state, void *dest, size_t len, Error *err);

extern int
initSeqReaderSet(SeqReaderSet *readerSet, int initialSuperSet,
                 int numConsumers, int *tags, SeqDataTranslator xltors[],
                 SeqDataReader *generatedReaders, size_t seqElemSize,
                 generatorFunc generator, void *generatorState)
{
  int i;
  initEmptySeqReaderSet(readerSet, initialSuperSet, seqElemSize,
                        generator, generatorState);
  for (i = 0; i < numConsumers; ++i)
  {
    generatedReaders[i]
      = seqReaderSetRegisterConsumer(readerSet, tags[i], xltors[i]);
    if (!generatedReaders[i].readData)
      break;
  }
  return i;
}

extern void
initEmptySeqReaderSet(SeqReaderSet *readerSet, int initialSuperSet,
                      size_t seqElemSize, generatorFunc generator,
                      void *generatorState)
{
  readerSet->tagSuperSet = initialSuperSet;
  readerSet->numConsumers = 0;
  readerSet->consumerList = NULL;
  readerSet->backlogStartPos = 0;
  readerSet->backlogLen = readerSet->backlogSize = 0;
  readerSet->backlogElemSize = seqElemSize;
  readerSet->seqDataBacklog = NULL; /* only allocate backlog when it's
                                   * really needed! */
  readerSet->generatorState = generatorState;
  readerSet->generator = generator;
}

extern void
destructSeqReaderSet(SeqReaderSet *readerSet)
{
  struct seqReaderState *p;
  assert(readerSet);
  ListDo(struct seqReaderState, readerSet->consumerList, ma_free(p));
  if (readerSet->autoConsumerList)
    ma_free(readerSet->autoConsumerList);
  if (readerSet->seqDataBacklog)
    ma_free(readerSet->seqDataBacklog);
}

extern SeqDataReader
seqReaderSetRegisterConsumer(SeqReaderSet *readerSet, int tag,
                             SeqDataTranslator xltor)
{
  int availId = readerSet->numConsumers++;
  SeqDataReader reader;
  struct seqReaderState *state;
  state = ma_malloc(sizeof (*state));
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

static size_t
seqReaderSetRead(void *src, void *dest, size_t len, Error *err)
{
  struct seqReaderState *state = src;
  struct seqReaderSet *readerSet;
  Seqpos pos;
  size_t elemsLeft;
  assert(src);
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
                             state->xltor, err);
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
      Seqpos subLen
        = MIN(elemsLeft, readerSet->backlogStartPos - pos
              + readerSet->backlogLen),
        charsWritten =
        SDRTranslate(state->xltor, dest, (char *)readerSet->seqDataBacklog
                     + (pos - readerSet->backlogStartPos)
                     * readerSet->backlogElemSize,
                     subLen);
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
static inline Seqpos
seqReaderSetFindMinOpenRequestBySet(SeqReaderSet *readerSet, int setSpec)
{
  Seqpos min = -1;
  int i;
  assert(readerSet);
  for (i = 0; i < readerSet->numConsumers; ++i)
    if (readerSet->consumers[i].tag & setSpec)
      min = MIN(min, readerSet->consumers[i].nextReadPos);
  return min;
}

static inline Seqpos
seqReaderSetFindMinOpenRequestByVal(SeqReaderSet *readerSet, int tagVal)
{
  Seqpos min = -1;
  int i;
  assert(readerSet);
  for (i = 0; i < readerSet->numConsumers; ++i)
    if (readerSet->consumers[i].tag == tagVal)
      min = MIN(min, readerSet->consumers[i].nextReadPos);
  return min;
}
#endif

static inline Seqpos
seqReaderSetFindMinOpenRequest(SeqReaderSet *readerSet)
{
  struct seqReaderState *p;
  Seqpos min = -1;
  assert(readerSet);
  p = readerSet->consumerList;
  while (p)
  {
    min = MIN(min, p->nextReadPos);
    p = p->next;
  }
  return min;
}

static void
seqReaderSetMove2Backlog(void *backlogState, const void *seqData,
                         Seqpos requestStart, size_t requestLen)
{
  Seqpos requestMinPos;
  struct seqReaderSet *readerSet = backlogState;
  assert(backlogState && (requestLen?(seqData!=NULL):1));
  requestMinPos = seqReaderSetFindMinOpenRequest(readerSet);
  /* 1. move still unread old values as far as possible to head of copy */
  assert(requestMinPos >= readerSet->backlogStartPos);
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
    Seqpos copyStartPos = MAX(requestMinPos, requestStart);
    size_t copyLen = requestLen - copyStartPos + requestStart;
    /* 2. extend backlog to also accept all values in last region read */
    if (copyLen)
    {
      size_t backlogSizeLeft
        = readerSet->backlogSize - readerSet->backlogLen;
      if (copyLen > backlogSizeLeft)
      {
        size_t newSize = readerSet->backlogLen + copyLen;
        readerSet->seqDataBacklog
          = ma_realloc(readerSet->seqDataBacklog, readerSet->backlogElemSize
                       * newSize);
        readerSet->backlogSize = newSize;
      }
      memcpy((char *)readerSet->seqDataBacklog
             + readerSet->backlogLen * readerSet->backlogElemSize,
             (char *)seqData
             + (copyStartPos - requestStart) * readerSet->backlogElemSize,
             readerSet->backlogElemSize * copyLen);
      readerSet->backlogLen += copyLen;
    }
  }
}

static inline int
seqReaderSetGetConsumerTag(struct seqReaderState *state)
{
  assert(state);
  return state->tag;
}

static inline void
seqReaderSetAdvanceConsumer(struct seqReaderState *state, size_t len)
{
  assert(state);
  state->nextReadPos += len;
}
