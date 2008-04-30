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

typedef struct suffixarrayFileInterface SuffixarrayFileInterface;

extern void
initSuffixarrayFileInterface(SuffixarrayFileInterface *sai,
                             Suffixarray *sa);

extern void
destructSuffixarrayFileInterface(SuffixarrayFileInterface *sai);

extern SeqDataReader
SAIMakeReader(SuffixarrayFileInterface *sai, enum sfxDataRequest rtype);

extern SeqDataReader
SAIMakeBWTReader(SuffixarrayFileInterface *sai);

extern SeqDataReader
SAIMakeSufTabReader(SuffixarrayFileInterface *sai);

/**
 * @brief Gets symbols of original sequence at given position.
 * @param state SuffixarrayFileInterface reference
 * @param dest write symbols here
 * @param pos get symbols starting at this position in original sequence
 * @param len length of string to read
 * @return actual number of symbols read
 */
extern size_t
SAIGetOrigSeqSym(void *state, Symbol *dest, Seqpos pos, size_t len);

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
 * \brief Get reference for original sequence object.
 *
 * @param sai SuffixarrayFileInterface reference
 * @return reference of sequence object
 */
static inline const Encodedsequence *
SAIGetEncSeq(const SuffixarrayFileInterface *sai);

/**
 * @brief Query appropriate alphabet encoding for suffix array.
 * @param state reference of Suffixarray object
 * @return alphabet
 */
extern MRAEnc *
newMRAEncFromSA(const Suffixarray *sa);

static inline MRAEnc *
newMRAEncFromSAI(const SuffixarrayFileInterface *sai);

/* visible for the compiler, but not meant for users to depend upon */
#include "eis-suffixarray-interface-priv.h"

#endif
