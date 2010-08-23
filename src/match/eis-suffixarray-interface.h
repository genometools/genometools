/*
  Copyright (C) 2007,2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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
 * Defines an interface deriving the abstract interface defined in
 * eis-sa-common.h for suffix array class objects (type SASeqSrc).
 */
#include <stdlib.h>
#include "core/error.h"
#include "match/sarr-def.h"

#include "match/eis-mrangealphabet.h"
#include "match/eis-seqdatasrc.h"
#include "match/eis-sequencemultiread.h"

typedef struct suffixarrayFileInterface SuffixarrayFileInterface;

void
gt_initSuffixarrayFileInterface(SuffixarrayFileInterface *sai,
                             unsigned long seqLen,
                             Suffixarray *sa);

SuffixarrayFileInterface *
gt_newSuffixarrayFileInterface(Suffixarray *sa, unsigned long seqLen);

void
gt_destructSuffixarrayFileInterface(SuffixarrayFileInterface *sai);

void
gt_deleteSuffixarrayFileInterface(SuffixarrayFileInterface *sai);

SeqDataReader
gt_SAIMakeReader(SuffixarrayFileInterface *sai, enum sfxDataRequest rtype);

SeqDataReader
gt_SAIMakeBWTReader(SuffixarrayFileInterface *sai);

SeqDataReader
gt_SAIMakeSufTabReader(SuffixarrayFileInterface *sai);

/**
 * @brief Gets symbols of original sequence at given position.
 * @param state SuffixarrayFileInterface reference
 * @param dest write symbols here
 * @param pos get symbols starting at this position in original sequence
 * @param len length of string to read
 * @return actual number of symbols read
 */
size_t
gt_SAIGetOrigSeq(const void *state, Symbol *dest, unsigned long pos,
                 size_t len);

/**
 * @brief Query position of suffix starting at position 0, can be
 * undefined if not yet encountered.
 *
 * @param state reference of Suffixarray object
 * @return
 */
Definedunsignedlong
gt_SAIGetRot0Pos(const void *state);

/**
 * \brief Get reference for original sequence object.
 *
 * @param sai SuffixarrayFileInterface reference
 * @return reference of sequence object
 */
static inline const GtEncseq *
SAIGetEncSeq(const SuffixarrayFileInterface *sai);

static inline GtReadmode
SAIGetGtReadmode(const SuffixarrayFileInterface *sai);

static inline unsigned long
SAIGetLength(const SuffixarrayFileInterface *sai);

/**
 * @brief Query appropriate alphabet encoding for suffix array.
 * @param state reference of Suffixarray object
 * @return alphabet
 */
MRAEnc *
gt_SANewMRAEnc(const GtAlphabet *gtalphabet);

static inline MRAEnc *
SAINewMRAEnc(const SuffixarrayFileInterface *sai);

static inline struct SASeqSrc *
SAI2SASS(SuffixarrayFileInterface *sai);

/* visible for the compiler, but not meant for users to depend upon */
#include "eis-suffixarray-interface-siop.h"

#endif
