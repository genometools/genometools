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

#ifndef EIS_SEQUENCEMULTIREAD_H
#define EIS_SEQUENCEMULTIREAD_H
/**
 * \file eis-sequencemultiread.h
 * Keeps information about multiple synchronous reads of sequence
 * data, where only one source and multiple consumers exist.
 */
#include "core/assert_api.h"
#include "core/minmax.h"
#include "match/eis-seqdatasrc.h"

/** every reader is identified by a unique scalar */
typedef unsigned consumerID;

/* moves data that can no longer be regenerated to the backlog */
typedef void (*move2BacklogFunc)(void *backlogState, const void *seqData,
                                 unsigned long requestStart, size_t requestLen);

/**
 * basic idea: let this function write the required data to output
 * and also call the move2Backlog callback if any data will be
 * invalidated after this call and is still required by other consumers
 * @return number of elements actually generated (might be short on eof etc.)
 */
typedef size_t (*generatorFunc)(void *generatorState, void *backlogState,
                                move2BacklogFunc move2Backlog, void *output,
                                unsigned long generateStart, size_t len,
                                SeqDataTranslator xltor);

typedef struct seqReaderSet SeqReaderSet;

struct seqReaderSet
{
  int numConsumers, numAutoConsumers;
  int tagSuperSet;
  struct seqReaderState *consumerList;
  struct seqSinkState *autoConsumerList;
  unsigned long backlogStartPos;
  bool fromSuffixsortspace;
  size_t backlogSize, backlogLen, backlogElemSize;
  void *seqDataBacklog, *generatorState;
  generatorFunc generator;
};

/**
 * @return numReaders if all consumers registered otherwise
 */
/*
int
gt_initSeqReaderSet(SeqReaderSet *readerSet, int initialSuperSet,
                 int numConsumers, int *tags, SeqDataTranslator xltors[],
                 SeqDataReader *generatedReaders, size_t seqElemSize,
                 generatorFunc generator, void *generatorState);
*/

void
gt_initEmptySeqReaderSet(SeqReaderSet *readerSet, int initialSuperSet,
                         bool fromSuffixsortspace,
                         size_t seqElemSize,
                         generatorFunc generator,
                         void *generatorState);

/**
 * @return readData field will be NULL on error -> test with
 * SDRIsValid */
SeqDataReader
gt_seqReaderSetRegisterConsumer(SeqReaderSet *readerSet, int tag,
                             SeqDataTranslator xltor);

/**
 * @brief The registered writer will be called automatically for any
 * data that is to be invalidated.
 *
 * @return false on error, true if successfully registered
 */
bool
gt_seqReaderSetRegisterAutoConsumer(SeqReaderSet *readerSet, int tag,
                                 SeqDataWriter writer);

void
gt_destructSeqReaderSet(SeqReaderSet *readerSet);

#endif
