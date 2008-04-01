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

#ifndef EIS_SUFFIXARRAY_INTERFACE_H
#define EIS_SUFFIXARRAY_INTERFACE_H

/**
 * \file eis-suffixarray-interface.h
 * Defines functions conforming to the signatures defined in
 * eis-construction-interface.h for suffix array objects.
 */
#include <stdlib.h>
#include "libgtcore/error.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-seqdatasrc.h"
#include "libgtmatch/eis-sequencemultiread.h"

/**
 * Used to pass suffixarray to read from and alphabet to encode with
 * to readers.
 */
struct suffixarrayFileInterface
{
  Suffixarray *sa;              /**< the suffix array to read from */
  struct seqReaderSet readerSet;
  int numBWTFileReaders;
  struct saTaggedXltorStateList xltorStates;
};

extern void
initSuffixarrayFileInterface(struct suffixarrayFileInterface *sai,
                             Suffixarray *sa);

extern void
destructSuffixarrayFileInterface(struct suffixarrayFileInterface *sai);

extern SeqDataReader
SAIMakeReader(struct suffixarrayFileInterface *sai, enum sfxDataRequest rtype);

extern SeqDataReader
SAIMakeBWTReader(struct suffixarrayFileInterface *sai);

extern SeqDataReader
SAIMakeSufTabReader(struct suffixarrayFileInterface *sai);

/**
 * @brief Gets symbols of original sequence at given position.
 * @param state reference of a struct suffixarrayFileInterface
 * @param dest write symbols here
 * @param pos get symbols starting at this position in original sequence
 * @param len length of string to read
 * @return actual number of symbols read
 */
extern size_t
SAIGetOrigSeqSym(void *state, Symbol *dest, Seqpos pos, size_t len);

#if 0
/**
 * @brief Read part of the suffix array starting from position after
 * last read.
 * @param src reference of a Suffixarray object
 * @param dest write suffix array values here
 * @param len length of part to read
 * @param err genometools error object reference
 * @return number of entries read, less than len if end of sequence
 * reached
 */
extern size_t
saReadSeqpos(void *src, Seqpos *dest, size_t len, Error *err);
#endif

/**
 * @brief Query position of suffix starting at position 0, can be
 * undefined if not yet encountered.
 *
 * @param state reference of Suffixarray object
 * @return
 */
extern DefinedSeqpos
reportSAILongest(void *state);

/**
 * @brief Query appropriate alphabet encoding for suffix array.
 * @param state reference of Suffixarray object
 * @return alphabet
 */
extern MRAEnc *
newMRAEncFromSA(const Suffixarray *sa);

static inline MRAEnc *
newMRAEncFromSAI(const struct suffixarrayFileInterface *sai)
{
  return newMRAEncFromSA(sai->sa);
}

#endif
