/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#ifndef EIS_BWTSEQCREATE_H
#define EIS_BWTSEQCREATE_H

#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseqpriv.h"
#include "libgtmatch/eis-encidxseq.h"

typedef int (*srcReadFunc)(Seqpos *dest, void *src, Env *env);
typedef struct encIdxSeq *(*indexCreateFunc)
  (void *src, Seqpos totalLen, const Str *projectName,
   const union bwtSeqParam *params, size_t numExtHeaders, uint16_t *headerIDs,
   uint32_t *extHeaderSizes, headerWriteFunc *extHeaderCallbacks,
   void **headerCBData, bitInsertFunc biFunc, BitOffset cwBitsPerPos,
   BitOffset maxBitsPerPos, void *cbState, Env *env);

extern EISeq *
createBWTSeqGeneric(const struct bwtParam *params, void *baseSrc,
                    indexCreateFunc createIndex, unsigned aggInterval,
                    srcReadFunc readCallback, Seqpos totalLen, Env *env);

#endif
