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
#include <assert.h>
#include "libgtcore/minmax.h"
#include "libgtmatch/eis-seqdatasrc.h"
#include "libgtmatch/seqpos-def.h"

/** every reader is identified by a unique scalar */
typedef unsigned consumerID;

/* moves data that can no longer be regenerated to the backlog */
typedef void (*move2BacklogFunc)(void *backlogState, const void *seqData,
                                 Seqpos requestStart, size_t requestLen);

/**
 * basic idea: let this function write the required data to output
 * and also call the move2Backlog callback if any data will be
 * invalidated after this call
 * @return number of elements actually generated (might be short on eof etc.)
 */
typedef size_t (*generatorFunc)(void *generatorState, void *backlogState,
                                move2BacklogFunc move2Backlog, void *output,
                                Seqpos generateStart, size_t len,
                                SeqDataTranslator xltor, Error *err);

typedef struct seqReaderSet SeqReaderSet;

struct seqReaderSet
{
  int numConsumers;
  int tagSuperSet;
  struct seqReaderState *consumerList;
  Seqpos backlogStartPos;
  size_t backlogSize, backlogLen, backlogElemSize;
  void *seqDataBacklog, *generatorState;
  generatorFunc generator;
};

/**
 * @return numReaders if all consumers registered otherwise
 */
extern int
initSeqReaderSet(SeqReaderSet *readerSet, int initialSuperSet,
                 int numConsumers, int *tags, SeqDataTranslator xltors[],
                 SeqDataReader *generatedReaders, size_t seqElemSize,
                 generatorFunc generator, void *generatorState);

extern void
initEmptySeqReaderSet(SeqReaderSet *readerSet, int initialSuperSet,
                      size_t seqElemSize, generatorFunc generator,
                      void *generatorState);

/**
 * @return readData field will be NULL on error */
extern SeqDataReader
seqReaderSetRegisterConsumer(SeqReaderSet *readerSet, int tag,
                             SeqDataTranslator xltor);

extern void
destructSeqReaderSet(SeqReaderSet *readerSet);

#endif
